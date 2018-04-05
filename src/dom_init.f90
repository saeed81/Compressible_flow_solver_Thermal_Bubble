MODULE dom_init
  USE phys_cst
  USE dcn_vctr
  USE dcn_opt
  USE step
  USE iom
  USE io_ezcdf
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: istate_init

CONTAINS 

  SUBROUTINE istate_init

    CALL test_namelist
    CALL get_cexper
    CALL get_config
    CALL comp_rgamma
    CALL read_lapsrt

    CALL exner_init
    CALL rho_init
    CALL alpha_init
    CALL perturb_init
    CALL par_step
    CALL name_file
    CALL write_init (nt)

  END SUBROUTINE istate_init

  SUBROUTINE exner_init

    INTEGER :: ji, jj

    DO jj = 1, nz
       TP(:,jj) = rn_lapsrate * (jj - 1) * dz + rthetas
    END DO

    EXNER(:,1) = rcp * ( rps / rps ) ** ( rd / rcp)

    IF( ln_neutral) THEN
       DO jj = 2, nz

          EXNER(:,jj) = EXNER(:,jj-1) - ( grav / TP(:,jj-1) ) * dz

       END DO
    ELSE
       DO jj = 2, nz

          EXNER(:,jj) = EXNER(:,jj-1) - ( grav * dz ) / (TP(:,jj) - &
               & TP(:,jj-1)) * LOG( TP(:, jj) / TP(:, jj-1))

       END DO
    END IF
    THETA(:,:) = TP(:,:)

  END SUBROUTINE exner_init

  SUBROUTINE rho_init

    P(:,1) = rps
    P(:,2:nz) = rps * ( (EXNER(:,2:nz) / rcp) ** ( rcp / rd ))
    PH(:,:) = P(:,:)
    RHO(:,:) = ( P(:,:) / ( rd * TP(:,:))) * ((rps / P(:,:)) ** (rd / rcp))

  END SUBROUTINE rho_init

  SUBROUTINE alpha_init
    INTEGER :: ji
    ALPHA(:) = 0.0_wp
    !
    DO ji = 1, nx
       GZ(:) = P(ji,:)
       !
       CALL forwardz( dz, GZ )
       !
       RF(:) = -GZ(:) / grav
       GZ(:) = P(ji,:)
       !
       CALL backwardz( dz, GZ )
       !
       RB(:) = - GZ(:) / grav
       ALPHA(2:nz) = ( RB(2:nz) - RHO(ji,1:nz-1) ) / ( RF(1:nz-1) + &
            & RB(2:nz) - 2.0_wp * RHO(ji,1:nz-1) )
    END DO

  END SUBROUTINE alpha_init

  SUBROUTINE perturb_init

    SELECT CASE (trim(cn_exper))

    CASE('straka')

       CALL straka_init

    CASE('wbrbb')

       CALL wbrbb_init

    CASE('wbr4b')

       CALL wbr4b_init

    CASE('dncrt')

       CALL dncrt_init   

    END SELECT

    W(:,:) = 0.0_wp
    U(:,:) = 0.0_wp

    RHOTHTA(:,:) = ( rps / rd ) * ((P(:,:) / rps) ** (rcv / rcp ))

    V(:,:,1) = RHO(:,:)
    V(:,:,2) = RHO(:,:) * U(:,:)
    V(:,:,3) = RHO(:,:) * W(:,:)
    V(:,:,4) = RHOTHTA(:,:)

  END SUBROUTINE perturb_init

  SUBROUTINE straka_init

    INTEGER :: ji, jj
    REAL(wp):: zvalue

    DO jj = 1, nz
       DO ji = 1, nx
          !
          CALL bstraka( (ji - 1) * dx, (jj - 1) * dz, zvalue)
          !
          zvalue = zvalue * rcp / EXNER(ji,jj)

          IF( zvalue <= 0.0_wp ) THEN
             THETA(ji,jj) = zvalue + TP(ji,jj)
             RHO(ji,jj) = ( P(ji,jj) / (rd * THETA(ji,jj))) * &
                  &(( rps / P(ji,jj)) ** ( rd / rcp ) )
          END IF

       END DO
    END DO

  END SUBROUTINE straka_init

  SUBROUTINE wbrbb_init
    INTEGER :: ji, jj
    REAL(wp):: zvalue


    DO jj = 1, nz
       DO ji = 1, nx
          CALL bwbrbb( (ji - 1) * dx, (jj - 1) * dz, zvalue)


          IF( zvalue >= 0.0_wp ) THEN
             PRINT*, zvalue
             THETA(ji,jj) = zvalue + TP(ji,jj)
             RHO(ji,jj) = ( P(ji,jj) / (rd * THETA(ji,jj))) * &
                  &(( rps / P(ji,jj)) ** ( rd / rcp ) )
          END IF

       END DO
    END DO
    
  END SUBROUTINE wbrbb_init

  SUBROUTINE wbr4b_init

    INTEGER :: ji, jj
    REAL(wp):: zvalue

    DO jj = 1, nz
       DO ji = 1, nx

          CALL bwbr4b( (ji - 1) * dx, (jj - 1) * dz, zvalue)

          IF( zvalue >= 0.0_wp ) THEN
             THETA(ji,jj) = zvalue + TP(ji,jj)
             RHO(ji,jj) = ( P(ji,jj) / (rd * THETA(ji,jj))) * &
                  &(( rps / P(ji,jj)) ** ( rd / rcp ) )
          END IF

       END DO
    END DO

  END SUBROUTINE wbr4b_init

  SUBROUTINE dncrt_init

    INTEGER :: ji, jj
    REAL(wp):: zvalue

    DO jj = 1, nz
       DO ji = 1, nx
          CALL bdncrt( (ji - 1) * dx, (jj - 1) * dz, zvalue)


          IF( zvalue <= 0.0_wp ) THEN
             THETA(ji,jj) = zvalue + TP(ji,jj)
             RHO(ji,jj) = ( P(ji,jj) / (rd * THETA(ji,jj))) * &
                  &(( rps / P(ji,jj)) ** ( rd / rcp ) )
          END IF

       END DO
    END DO

  END SUBROUTINE dncrt_init

  SUBROUTINE bstraka( px, pz, pvalue )

    REAL(wp), INTENT(in  ) :: px, pz 
    REAL(wp), INTENT( out) :: pvalue
    REAL(wp)               :: xr, zr
    REAL(wp)               :: xc, zc
    REAL(wp)               :: arx, arz
    REAL(wp)               :: beta

    xc = 0.0_wp
    zc = 3000.0_wp
    xr = 4000.0_wp
    zr = 2000.0_wp

    arx = ( px - xc ) / xr
    arz = ( pz - zc ) / zr
    beta = SQRT ( arx * arx + arz * arz )

    IF( beta <= 1.0_wp )  THEN
       pvalue = -15.0_wp * ( COS ( rpi * beta ) + 1.0_wp ) / 2.0_wp
    ELSE
       pvalue = 1.0_wp
    END IF

  END SUBROUTINE bstraka

  SUBROUTINE bwbr4b(px, pz, pvalue)

    REAL(wp), INTENT(in  ) :: px, pz 
    REAL(wp), INTENT( out) :: pvalue
    REAL(wp)               :: xc, zc
    REAL(wp)               :: rb, rdst
    xc = 1.6e3_wp
    zc = 1.e3_wp
    rb = 1.0e3_wp

    rdst = SQRT((px - xc)* (px -xc) + &
         & (pz - zc) * (pz - zc))
    IF( rdst <= rb ) THEN
       pvalue = -2.0_wp * ( rdst / rb ) + 2.0_wp
    ELSE
       pvalue = -1.0_wp
    END IF

  END  SUBROUTINE bwbr4b

  SUBROUTINE bwbrbb ( px, pz, pvalue )

    REAL(wp), INTENT(in  ) :: px, pz 
    REAL(wp), INTENT( out) :: pvalue
    REAL(wp)               :: xr, zr
    REAL(wp)               :: xc, zc
    REAL(wp)               :: arx, arz
    REAL(wp)               :: beta

    xc = 20000.0_wp
    zc = 2750.0_wp
    xr = 2500.0_wp
    zr = xr

    arx = ( px - xc ) / xr 
    arz = ( pz - zc ) / zr

    beta = SQRT ( arx * arx + arz * arz )

    IF( beta <= 1.0_wp )  THEN

       pvalue = 6.60_wp * COS( rpi * beta / 2.0_wp) &
            & * COS ( rpi * beta / 2.0_wp ) 
    ELSE
       pvalue = -1.0_wp
    END IF

  END SUBROUTINE bwbrbb


  SUBROUTINE bdncrt(px, pz, pvalue)

    REAL(wp), INTENT(in  ) :: px, pz 
    REAL(wp), INTENT( out) :: pvalue
    REAL(wp)               :: xc, zc

    xc = 5000.0_wp
    zc = 5000.0_wp
    !
    IF((ABS(px) <= xc) .AND. (pz <= zc)) THEN
       pvalue = 10.0_wp * ( pz / zc )- 10.0_wp
    ELSE
       pvalue = 2.0_wp
    END IF

  END SUBROUTINE bdncrt

  SUBROUTINE test_namelist

    LOGICAL               ::   llok      ! check the existence

    INQUIRE( FILE = "namelist", EXIST = llok )

    IF( .NOT.llok ) THEN
       STOP "please provide the namelist first"
    END IF
    OPEN( unit=numnam, FILE="namelist", FORM="FORMATTED", STATUS="OLD" )

  END SUBROUTINE test_namelist

  SUBROUTINE write_init (kt)
    INTEGER :: ji
    INTEGER, INTENT(in) :: kt
    
    ALLOCATE(vtime(kt))

    DO ji = 1, nx
       vlon(ji) = (ji - 1) * dx
    END DO
    
    DO ji = 1, nz
       vht(ji) = (ji - 1) * dz
    END DO
    
       DO ji = 1, kt
          vtime(ji) = (ji - 1) * rn_twrite
       END DO


    CALL P2D_T(idx_f1, idx_v1, nx, nz, kt, 1, real(vlon,kind=4), real(vht,kind=4),&
         & real(vtime,kind=4), real(U,kind=4), trim(cfilu), &
         & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_u), trim(cunit_u), &
         & trim(cln_u), 0.0)

    CALL P2D_T(idx_f2, idx_v2, nx, nz, kt, 1, real(vlon,kind=4), real(vht,kind=4),&
         & real(vtime,kind=4), real(W,kind=4), trim(cfilw), &
         & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_w), trim(cunit_w), &
         & trim(cln_w), 0.0)

    CALL P2D_T(idx_f3, idx_v3, nx, nz, kt, 1, real(vlon,kind=4), real(vht,kind=4),&
         & real(vtime,kind=4), real(THETA-TP,kind=4), trim(cfilpt), &
         & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_pt), trim(cunit_pt), &
         & trim(cln_pt), 0.0)

    CALL P2D_T(idx_f4, idx_v4, nx, nz, kt, 1, real(vlon,kind=4), real(vht,kind=4),&
         & real(vtime,kind=4), real(P-PH,kind=4), trim(cfilp), &
         & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_p), trim(cunit_p), &
         & trim(cln_p), 0.0)


    CALL P2D_T(idx_f5, idx_v5, nx, nz, kt, 1, real(vlon,kind=4), real(vht,kind=4),&
         & real(vtime,kind=4), real(ZETA,kind=4), trim(cfilvor), &
         & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_vor), trim(cunit_vor), &
         & trim(cln_vor), 0.0)

    CALL P2D_T(idx_f6, idx_v6, nx, nz, kt, 1, real(vlon,kind=4), real(vht,kind=4),&
         & real(vtime,kind=4), real(DIV,kind=4), trim(cfildiv), &
         & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_div), trim(cunit_div), &
         & trim(cln_div), 0.0)


  END SUBROUTINE write_init

END MODULE dom_init
