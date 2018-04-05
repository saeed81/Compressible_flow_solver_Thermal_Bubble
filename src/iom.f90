MODULE iom
  !!====================================================================================
  !!                       ***  MODULE  iom  ***
  !!               Declaration of logical units for all the files 
  !!                used as either in or out in dcn model
  !!====================================================================================
  USE par_dcn
  USE in_out_manager
  USE step

  IMPLICIT NONE
  PUBLIC

  PUBLIC :: name_file


  !!-------------------------------------------------------------------------------------


CONTAINS

  SUBROUTINE name_file

    INTEGER  ::   immss, inbmn, indg
    INTEGER :: jg, ji
    REAL(wp)  ::  rmmss = 60.0
    CHARACTER(len=100) :: clfmt0, clfmt, clave
    CHARACTER(len=100)::  clrun, clftmp
    CHARACTER(len=3), DIMENSION(6)   :: cldgrid
    CHARACTER(len=4), DIMENSION(2)   :: clresn
    CHARACTER(len=4)                 :: clrtmp


    NAMELIST/namcdout/cn_drout

    REWIND(numnam)
    READ(numnam,namcdout)

    cldgrid = (/"U  ","W  ","T  ","P  ","VOR","DIV"/)



    WRITE(clrtmp,'(I4.4)') nx
    IF(clrtmp(1:1) == "0" ) THEN
       IF(clrtmp(2:2) == "0" ) THEN
          WRITE(clresn(1),'(I2.2)')nx
       ELSE
          WRITE(clresn(1),'(I3.3)')nx
       END IF
    ELSE
       WRITE(clresn(1),'(I4.4)')nx 
    END IF

    WRITE(clrtmp,'(I4.4)')nz
    IF(clrtmp(1:1) == "0" ) THEN
       IF(clrtmp(2:2) == "0" ) THEN
          WRITE(clresn(2),'(I2.2)')nz
       ELSE
          WRITE(clresn(2),'(I3.3)')nz
       END IF
    ELSE
       WRITE(clresn(2),'(I4.4)')nz 
    END IF



    IF(len_trim(cn_drout) /= index(cn_drout,"/",back = .TRUE.))THEN
       cn_drout=TRIM(cn_drout)//"/"

    END IF


    DO jg = 1, SIZE(cldgrid)
       immss = NINT( rmmss )
       clfmt0 = "('(a,i',i1,',a)')"

       IF( MOD( nn_twrite, immss ) == 0 ) THEN     ! frequency in minutes
          inbmn  = nn_twrite / immss
          indg   = INT(LOG10(REAL(inbmn ,wp))) + 1  ! number of digits needed to write minutes freque\

          WRITE(clfmt, clfmt0) indg

          WRITE(clave, clfmt) '_', inbmn , 'mn'
       ELSE
          indg   = INT(LOG10(REAL(nn_twrite,wp))) + 1   ! number of digits needed to write seconds frequency
          WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', nn_twrite, 's'
       ENDIF



       immss = NINT( rmmss )
       clfmt0 = "('(a,i',i1,',a)')"

       IF( MOD( nn_twrite, immss ) == 0 ) THEN


          inbmn  = nn_tmax / immss
          indg   = INT(LOG10(REAL(inbmn ,wp))) + 1  ! number of digits needed to write minutes freque\

          WRITE(clfmt, clfmt0) indg

          WRITE(clrun, clfmt) '_', inbmn , 'mn'
       ELSE
          indg   = INT(LOG10(REAL(nn_tmax,wp))) + 1   ! number of digits needed to write seconds frequency

          WRITE(clfmt, clfmt0) indg                ;   WRITE(clrun, clfmt) '_', nn_tmax, 's'

       ENDIF
       
       nt = INT ( nn_tmax / nn_twrite) + 1
       

       WRITE(clftmp, '(a,a,"_",a,a,a,"_",a,"-",a,".nc")')TRIM(cn_drout), TRIM(cldgrid(jg)), TRIM(cn_exper),&
            &TRIM(clave),TRIM(clrun),TRIM(ADJUSTL(clresn(1))),TRIM(ADJUSTL(clresn(2)))

       SELECT CASE ( jg )
       CASE ( 1 )
          cfilu   = TRIM(clftmp)
          cvarout_u = "horizontal_velocity"
          cunit_u   = "ms-1"
          cln_u     = trim(cvarout_u)//"_"//TRIM(cn_exper)//"_"//TRIM(clave)//&
               TRIM(clrun)//TRIM(ADJUSTL(clresn(1)))//"-"//TRIM(ADJUSTL(clresn(2)))
       CASE ( 2 )
          cfilw   = TRIM(clftmp)
          cvarout_w = "vertical_velocity"
          cunit_w   = "ms-1"
          cln_w    = trim(cvarout_w)//"_"//TRIM(cn_exper)//"_"//TRIM(clave)//&
               TRIM(clrun)//TRIM(ADJUSTL(clresn(1)))//"-"//TRIM(ADJUSTL(clresn(2)))

       CASE ( 3 )
          cfilpt   = TRIM(clftmp)
          cvarout_pt  = "potential_temperature"
          cunit_pt   = "K"
          cln_pt    = trim(cvarout_pt)//"_"//TRIM(cn_exper)//"_"//TRIM(clave)//&
               TRIM(clrun)//TRIM(ADJUSTL(clresn(1)))//"-"//TRIM(ADJUSTL(clresn(2)))

       CASE ( 4 )
          cfilp  = TRIM(clftmp)
          cvarout_p = "pressure"
          cunit_p   = "Pa"
          cln_p     = trim(cvarout_p)//"_"//TRIM(cn_exper)//"_"//TRIM(clave)//&
               TRIM(clrun)//TRIM(ADJUSTL(clresn(1)))//"-"//TRIM(ADJUSTL(clresn(2)))

       CASE ( 5 )
          cfilvor = TRIM(clftmp)
          cvarout_vor = "vorticity"
          cunit_vor  = "s-1"
          cln_vor    = trim(cvarout_vor)//"_"//TRIM(cn_exper)//"_"//TRIM(clave)//&
               TRIM(clrun)//TRIM(ADJUSTL(clresn(1)))//"-"//TRIM(ADJUSTL(clresn(2)))

       CASE ( 6 )
          cfildiv = TRIM (clftmp)
          cvarout_div = "divergence"
          cunit_div   = "s-1"
          cln_div     = trim(cvarout_div)//"_"//TRIM(cn_exper)//"_"//TRIM(clave)//&
               TRIM(clrun)//TRIM(ADJUSTL(clresn(1)))//"-"//TRIM(ADJUSTL(clresn(2)))

       END SELECT

    END DO

  END SUBROUTINE name_file

END MODULE iom
