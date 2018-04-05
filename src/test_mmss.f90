PROGRAM test_mmss
  IMPLICIT NONE
  REAL(8)  :: rmmss=60.0
  INTEGER  :: immss, inbsec=240,inbmn,indg, inrsec=960
  CHARACTER(len=100)::clfmt0,clfmt,clave, clrun, cfile, cexper, cdrout
  CHARACTER(len=3), DIMENSION(6)   :: cldgrid
  CHARACTER(len=6),DIMENSION(4)    :: clexper
  CHARACTER(len=4), DIMENSION(2)   :: cresln
  CHARACTER(len=4)                 :: ctmp
  INTEGER :: jg, je, ji, nx, nz
  INTEGER, DIMENSION(2)            :: ifreq_run

  cldgrid = (/"U  ","W  ","T  ","P  ","VOR","DIV"/)
  clexper = (/"straka","wbrbb ","wbr4b ","dnctr "/)
  ifreq_run = (/inbsec, inrsec/)
  cdrout="/data/dcn/"
  cexper="straka"

  nx = 101
  nz = 301
  !!PRINT*, index(cdrout,"/",back = .TRUE.)
  !!PRINT*, len_trim(cdrout)
  WRITE(ctmp,'(I4.4)')nx
  IF(ctmp(1:1) == "0" ) THEN
     IF(ctmp(2:2) == "0" ) THEN
        WRITE(cresln(1),'(I2.2)')nx
     ELSE
        WRITE(cresln(1),'(I3.3)')nx
     END IF
  ELSE
     WRITE(cresln(1),'(I4.4)')nx 
  END IF

  WRITE(ctmp,'(I4.4)')nz
  IF(ctmp(1:1) == "0" ) THEN
     IF(ctmp(2:2) == "0" ) THEN
        WRITE(cresln(2),'(I2.2)')nz
     ELSE
        WRITE(cresln(2),'(I3.3)')nz
     END IF
  ELSE
     WRITE(cresln(2),'(I4.4)')nz 
  END IF



  IF(len_trim(cdrout) /= index(cdrout,"/",back = .TRUE.))THEN
     cdrout=trim(cdrout)//"/"
     !!PRINT*, cdrout
  END IF

  !!WRITE(*,*) 'nb of seconds per minute rmmss = ', rmmss, ' s'


  DO je = 1, SIZE(clexper)
     DO jg = 1, SIZE(cldgrid)
        immss = NINT( rmmss )
        clfmt0 = "('(a,i',i1,',a)')"

        IF( MOD( inbsec, immss ) == 0 ) THEN                              ! frequency in minutes
           inbmn  = inbsec / immss
           indg   = INT(LOG10(REAL(inbmn ,8))) + 1  ! number of digits needed to write minutes freque\

           WRITE(clfmt, clfmt0) indg
           !!print*,indg
           !!PRINT*,clfmt
           WRITE(clave, clfmt) '_', inbmn , 'mn'
        ELSE
           indg   = INT(LOG10(REAL(inbsec,8))) + 1   ! number of digits needed to write seconds frequency
           WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', inbsec, 's'
        ENDIF

        !!    PRINT*,trim(clave)

        immss = NINT( rmmss )
        clfmt0 = "('(a,i',i1,',a)')"

        IF( MOD( inbsec, immss ) == 0 ) THEN
           ! frequency in minutes
           !!IF( MOD( inrsec, immss ) == 0 ) THEN

           inbmn  = inrsec / immss
           indg   = INT(LOG10(REAL(inbmn ,8))) + 1  ! number of digits needed to write minutes freque\

           WRITE(clfmt, clfmt0) indg
           !!print*,indg
           !!PRINT*,clfmt
           WRITE(clrun, clfmt) '_', inbmn , 'mn'
        ELSE
           indg   = INT(LOG10(REAL(inrsec,8))) + 1   ! number of digits needed to write seconds frequency

           WRITE(clfmt, clfmt0) indg                ;   WRITE(clrun, clfmt) '_', inrsec, 's'

        ENDIF

        WRITE(cfile,'(a,a,"_",a,a,a,"_",a,"-",a,".PLT")')trim(cdrout),trim(cldgrid(jg)),trim(clexper(je)),&
             &trim(clave),trim(clrun),trim(adjustl(cresln(1))),trim(adjustl(cresln(2)))

        PRINT*,trim(cfile)
     END DO
  END DO

END PROGRAM test_mmss
