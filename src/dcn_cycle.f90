MODULE dcn_cycle
  USE phys_cst
  USE dcn_vctr
  USE iom
  USE dcn_opt
  USE step
  USE io_ezcdf
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cycl_mck

CONTAINS

  SUBROUTINE cycl_mck ( kt, ptcounter, kcycle )  

    INTEGER, INTENT(in)     :: kcycle
    INTEGER, INTENT(inout)  :: kt
    REAL(wp), INTENT(inout) :: ptcounter
    INTEGER :: ji, jj, jk, icount = 1



    E(:,:,1)   = U(:,:) * RHO(:,:)
    E(:,:,2)   = U(:,:) * RHO(:,:) * U(:,:) + P(:,:)
    E(:,:,3)   = U(:,:) * RHO(:,:) * W(:,:)
    E(:,:,4)   = U(:,:) * RHO(:,:) * THETA(:,:)
    F(:,:,1)   = W(:,:) * RHO(:,:)
    F(:,:,2)   = W(:,:) * RHO(:,:) * U(:,:)
    F(:,:,3)   = W(:,:) * RHO(:,:) * W(:,:)
    F(:,:,4)   = W(:,:) * RHO(:,:) * THETA(:,:)
    H(:,:,1)   = 0.0_wp
    H(:,:,2)   = 0.0_wp
    H(:,:,3)   = - RHO(:,:) * grav
    H(:,:,4)   = 0.0_wp
    PP(:,:,1)  = 0.0_wp
    PP(:,:,2)  = 0.0_wp
    PP(:,:,3)  = P(:,:)
    PP(:,:,4)  = 0.0_wp
    DIS(:,:,:) = 0.0_wp
    DPP(:,:,:) = 0.0_wp

    DO jk = 1, 4
       DO jj = 1, nz  

          GX(:) = E(:,jj,jk)

          IF ( kcycle == 1 .OR. kcycle == 3 ) THEN

             CALL forwardx ( dx, GX )
          ELSEIF ( kcycle == 2 .OR. kcycle == 4) THEN   
             CALL backwardx ( dx, GX )
          END IF


          DE(:,jj,jk) = GX(:)

       END DO
    END DO

    DO jk = 1, 4  
       DO ji = 1, nx 

          GZ(:) = F(ji,:,jk)

          IF ( kcycle == 1 .OR. kcycle == 4 ) THEN

             CALL backwardx ( dz, GZ )
          ELSE IF ( kcycle == 2 .OR. kcycle == 3) THEN
             CALL forwardx( dz, GZ ) 
          ENDIF


          DF(ji,:,jk) = GZ(:)

       END DO
    END DO
    DO ji = 1, nx 

       GZ(:) = PP(ji,:,3)

       IF ( kcycle == 1 .OR. kcycle == 4 ) THEN

          CALL backwardz ( dz, GZ )
       ELSE IF ( kcycle == 2 .OR. kcycle == 3 ) THEN
          CALL forwardz ( dz, GZ )
       END IF


       DPP(ji,:,3) = GZ(:)

    END DO

    CALL dissipation (dx, dz, U, RHO, rn_nu, D )


    DIS(2:nx-1,2:nz-1,2) = D(2:nx-1,2:nz-1)


    CALL dissipation( dx, dz, W, RHO, rn_nu, D)


    DIS(2:nx-1,2:nz-1,3) = D(2:nx-1,2:nz-1)

    CALL dissipation( dx, dz, THETA, RHO, rn_nu, D )

    DIS(2:nx-1,2:nz-1,4) = D(2:nx-1,2:nz-1)

    !********************************************************************  

    IF ( kcycle == 1 .OR. kcycle == 4 ) THEN

       Q(:,1,:) = V(:,1,:) - (rn_dt) * (DE(:,1,:)) - rn_dt * DF(:,1,:)

       Q(:,nz,:) = V(:,nz,:) - (rn_dt) * (DE(:,nz,:)) - (rn_dt) * (DF(:,nz,:)) + (rn_dt) * &
            & ((1._wp - alpha(nz)) * H(:,nz-1,:) + alpha(nz) * H(:,nz,:)) - &
            (rn_dt)*( DPP(:,nz,:) )

    ELSE IF ( kcycle == 2 .OR. kcycle == 3) THEN

       Q(:,1,:) = V(:,1,:) - (rn_dt) * (DE(:,1,:)) - (rn_dt) * (DF(:,1,:)) + (rn_dt) * &
            & ((1._wp - alpha(2)) * H(:,1,:) + alpha(2) * H(:,2,:)) - &
            (rn_dt)*( DPP(:,1,:) )          
       Q(:,nz,:) = V(:,nz,:) - (rn_dt) * (DE(:,nz,:)) - (rn_dt) * (DF(:,nz,:)) 

    END IF

    IF ( kcycle == 1 .OR. kcycle == 4 ) THEN

       DO jj = 2, nz-1
          Q(:,jj,:) = V(:,jj,:) - (rn_dt) * (DE(:,jj,:)) - (rn_dt) * (DF(:,jj,:)) + (rn_dt) &
               & * (alpha(jj) * H(:,jj,:) + ( 1.0_wp - alpha(jj) ) * &
               & H(:,jj-1,:)) - (rn_dt) * ( DPP(:,jj,:)) + rn_dt * DIS(:,jj,:)
       END DO

    ELSE IF ( kcycle == 2 .OR. kcycle == 3) THEN

       DO jj = 2, nz-1
          Q(:,jj,:) = V(:,jj,:) - (rn_dt) * (DE(:,jj,:)) - (rn_dt) * (DF(:,jj,:)) + (rn_dt) &
               & * (alpha(jj+1) * H(:,jj,:) + ( 1.0_wp - alpha(jj+1) ) * &
               & H(:,jj+1,:)) - (rn_dt) * ( DPP(:,jj,:)) + rn_dt * DIS(:,jj,:)
       END DO
    END IF



    IF (ln_rb) THEN
       Q(:,1,:) = 2.0_wp * Q(:,2,:) - Q(:,3,:)
    END IF

    IF (ln_ru) THEN
       Q(:,nz,:) = 2.0_wp * Q(:,nz-1,:) - Q(:,nz-2,:)
    END IF

    IF (ln_rl) THEN
       Q(1,:,:) = 2.0_wp * Q(2,:,:) - Q(3,:,:)
    END IF

    IF (ln_rr) THEN
       Q(nx,:,:) = 2.0_wp * Q(nx-1,:,:) - Q(nx-2,:,:)
    END IF

    RHOIN(:,:)   = Q(:,:,1)
    UIN(:,:)     = Q(:,:,2) / Q(:,:,1)
    WIN(:,:)     = Q(:,:,3) / Q(:,:,1)
    THETAIN(:,:) = Q(:,:,4) / Q(:,:,1)


    IF (ln_rb) THEN
       WIN(:,1) = 0.0_wp
    END IF
    IF ( ln_ru) THEN
       WIN(:,nz) = 0.0_wp
    END IF

    IF (ln_rl) THEN
       UIN(1,:) = 0.0_wp
    END IF
    IF (ln_rr) THEN
       UIN(nx,:) = 0.0_wp
    END IF

    PIN(:,:) = rps*( ( RHOIN(:,:) * THETAIN(:,:) * rd / rps ) ** rgamma )

    EIN(:,:,1)  = RHOIN(:,:) * UIN(:,:)
    EIN(:,:,2)  = UIN(:,:) * RHOIN(:,:) * UIN(:,:) + PIN(:,:)
    EIN(:,:,3)  = UIN(:,:) * RHOIN(:,:) * WIN(:,:)
    EIN(:,:,4)  = UIN(:,:) * RHOIN(:,:) * THETAIN(:,:)
    FIN(:,:,1)  = RHOIN(:,:) * WIN(:,:)
    FIN(:,:,2)  = WIN(:,:) * RHOIN(:,:) * UIN(:,:)
    FIN(:,:,3)  = WIN(:,:) * RHOIN(:,:) * WIN(:,:)
    FIN(:,:,4)  = WIN(:,:) * RHOIN(:,:) * THETAIN(:,:)
    HIN(:,:,1)  = 0.0_wp
    HIN(:,:,2)  = 0.0_wp
    HIN(:,:,3)  = - RHOIN(:,:) * grav
    HIN(:,:,4)  = 0.0_wp
    PPIN(:,:,1) = 0.0_wp
    PPIN(:,:,2) = 0.0_wp
    PPIN(:,:,3) = PIN(:,:)
    PPIN(:,:,4) = 0.0_wp
    DIS(:,:,:)  = 0.0_wp

    DO jk = 1, 4
       DO jj = 1, nz  

          GX(:) = EIN(:,jj,jk)

          IF( kcycle == 1 .OR. kcycle == 3) THEN
             CALL backwardx ( dx, GX ) 
          ELSE IF (kcycle == 2 .OR. kcycle == 4) THEN
             CALL forwardx ( dx, GX )
          END IF


          DE(:,jj,jk) = GX(:)

       END DO
    END DO
    DO jk = 1, 4  
       DO ji = 1, nx 

          GZ(:) = FIN(ji,:,jk)

          IF( kcycle == 1 .OR. kcycle == 4 ) THEN 
             CALL forwardx ( dz, GZ )
          ELSE IF( kcycle == 2 .OR. kcycle == 3 ) THEN 
             CALL backwardx ( dz, GZ )
          END IF

          DF(ji,:,jk) = GZ(:)

       END DO
    END DO

    DO ji = 1, nx 
       GZ(:) = PPIN(ji,:,3)

       IF( kcycle == 1 .OR. kcycle == 4 ) THEN
          CALL forwardz(dz,GZ)
       ELSE IF ( kcycle == 2 .OR. kcycle == 3 ) THEN
          CALL backwardz(dz,GZ)
       END IF

       DPP(ji,:,3) = GZ(:)
    END DO


    CALL dissipation( dx, dz, UIN, RHOIN, rn_nu, D )

    DIS(2:nx-1,2:nz-1,2) = D(2:nx-1,2:nz-1)

    CALL dissipation (dx, dz, WIN, RHOIN, rn_nu, D ) 

    DIS(2:nx-1,2:nz-1,3) = D(2:nx-1,2:nz-1)

    CALL dissipation( dx, dz, THETAIN, RHOIN, rn_nu, D ) 

    DIS(2:nx-1,2:nz-1,4) = D(2:nx-1,2:nz-1)



    IF(kcycle == 1 .OR. kcycle == 4) THEN 
       V(:,1,:) = (0.5_wp) * (V(:,1,:) + Q(:,1,:) - (rn_dt) * (DE(:,1,:))&
            & - (rn_dt) * (DF(:,1,:)) + (rn_dt) * ((1._wp- alpha(2)) * HIN(:,2,:) &
            & + alpha(2) * HIN(:,1,:)) - (rn_dt) * (DPP(:,1,:)))

       V(:,nz,:) = (0.5_wp) * (V(:,nz,:) + Q(:,nz,:) - (rn_dt) * (DE(:,nz,:)) &
            - (rn_dt) * (DF(:,nz,:)))
    END IF
    IF(kcycle == 2 .OR. kcycle == 3) THEN 

       V(:,1,:) = (0.5_wp) * (V(:,1,:) + Q(:,1,:) - (rn_dt) * (DE(:,1,:)) &
            - (rn_dt) * (DF(:,1,:)))  
       V(:,nz,:) = (0.5_wp) * (V(:,nz,:) + Q(:,nz,:) - (rn_dt) * (DE(:,nz,:)) &
            - (rn_dt) * (DF(:,nz,:))+ (rn_dt) * ((1.0_wp - alpha(nz)) * &
            & HIN(:,nz-1,:) + alpha(nz) * HIN(:,nz,:)) - (rn_dt) * (DPP(:,nz,:)))
    END IF

    IF(kcycle == 1 .OR. kcycle == 4) THEN 
       DO jj = 2, nz-1

          V(:,jj,:) = (0.5_wp) * (V(:,jj,:) + Q(:,jj,:) - (rn_dt) * (DE(:,jj,:))&
               & - (rn_dt) * (DF(:,jj,:)) + (rn_dt) * ((1.0_wp - alpha(jj+1)) * HIN(:,jj+1,:) + &
               & alpha(jj+1) * HIN(:,jj,:)) - (rn_dt) * (DPP(:,jj,:)) + rn_dt * DIS(:,jj,:))
       END DO
    END IF
    IF(kcycle == 2 .OR. kcycle == 3) THEN 
       DO jj = 2, nz-1

          V(:,jj,:) = (0.5_wp) * (V(:,jj,:) + Q(:,jj,:) - (rn_dt) * (DE(:,jj,:))&
               & - (rn_dt) * (DF(:,jj,:)) + (rn_dt) * ((1.0_wp - alpha(jj)) * HIN(:,jj-1,:) + &
               & alpha(jj) * HIN(:,jj,:)) - (rn_dt) * (DPP(:,jj,:)) + rn_dt * DIS(:,jj,:))
       END DO
    END IF

    IF (ln_rb) THEN
       V(:,1,:) = 2.0_wp * V(:,2,:) - V(:,3,:)
    END IF
    IF (ln_ru) THEN
       V(:,nz,:) = 2.0_wp * V(:,nz-1,:) -V(:,nz-2,:)
    END IF
    IF (ln_rl) THEN      
       V(1,:,:) = 2.0_wp * V(2,:,:) - V(3,:,:)
    END IF
    IF (ln_rr) THEN
       V(nx,:,:) = 2.0_wp * V(nx-1,:,:) - V(nx-2,:,:)
    END IF
    !
    kt = kt + 1
    ptcounter = kt * rn_dt
    !
    WRITE(numstp,*)ptcounter
    REWIND(numstp)
    !
    WRITE(numscr,*)ptcounter
    !
    RHO(:,:)   = V(:,:,1)
    U(:,:)     = V(:,:,2) / V(:,:,1)
    W(:,:)     = V(:,:,3) / V(:,:,1)
    THETA(:,:) = V(:,:,4) / V(:,:,1)
    !
    IF (ln_rb) THEN
       W(:,1) = 0.0_wp
    END IF
    IF ( ln_ru) THEN
       W(:,nz) = 0.0_wp
    END IF
    !
    IF (ln_rl) THEN
       U(1,:) = 0.0_wp
    END IF
    IF (ln_rr) THEN
       U(nx,:) = 0.0_wp
    END IF
    !
    P(:,:) = rps * ( ( RHO(:,:) * THETA(:,:) * rd / rps) ** rgamma )
    !
    E(:,:,1)   = U(:,:) * RHO(:,:)
    E(:,:,2)   = U(:,:) * RHO(:,:) * U(:,:) + P(:,:)
    E(:,:,3)   = U(:,:) * RHO(:,:) * W(:,:)
    E(:,:,4)   = U(:,:) * RHO(:,:) * THETA(:,:)
    F(:,:,1)   = W(:,:) * RHO(:,:)
    F(:,:,2)   = W(:,:) * RHO(:,:) * U(:,:)
    F(:,:,3)   = W(:,:) * RHO(:,:) * W(:,:)
    F(:,:,4)   = W(:,:) * RHO(:,:) * THETA(:,:)
    H(:,:,1)   = 0.0_wp
    H(:,:,2)   = 0.0_wp
    H(:,:,3)   = - RHO(:,:) * grav
    H(:,:,4)   = 0.0_wp
    PP(:,:,1)  = 0.0_wp
    PP(:,:,2)  = 0.0_wp
    PP(:,:,3)  = P(:,:)
    PP(:,:,4)  = 0.0_wp
    DIS(:,:,:) = 0.0_wp
    DPP(:,:,:) = 0.0_wp

    IF ( MOD ( ptcounter, rn_twrite ) == 0.0_wp) THEN

       CALL VORTICITY  ( dx, dz, U, W, ZETA )
       CALL DIVERGENCE ( dx, dz, U, W, DIV )
       icount = icount + 1
       PRINT*, 'icount', icount
       CALL P2D_T(idx_f1, idx_v1, nx, nz, nt, icount, real(vlon,kind=4), real(vht,kind=4),&
            & real(vtime,kind=4), real(U,kind=4), trim(cfilu), &
            & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_u), trim(cunit_u), &
            & trim(cln_u), 0.0)

       CALL P2D_T(idx_f2, idx_v2, nx, nz, nt, icount, real(vlon,kind=4), real(vht,kind=4),&
            & real(vtime,kind=4), real(W,kind=4), trim(cfilw), &
            & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_w), trim(cunit_w), &
            & trim(cln_w), 0.0)

       CALL P2D_T(idx_f3, idx_v3, nx, nz, nt, icount, real(vlon,kind=4), real(vht,kind=4),&
            & real(vtime,kind=4), real(THETA-TP,kind=4), trim(cfilpt), &
            & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_pt), trim(cunit_pt), &
            & trim(cln_pt), 0.0)

       CALL P2D_T(idx_f4, idx_v4, nx, nz, nt, icount, real(vlon,kind=4), real(vht,kind=4),&
            & real(vtime,kind=4), real(P-PH,kind=4), trim(cfilp), &
            & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_p), trim(cunit_p), &
            & trim(cln_p), 0.0)


       CALL P2D_T(idx_f5, idx_v5, nx, nz, nt, icount, real(vlon,kind=4), real(vht,kind=4),&
            & real(vtime,kind=4), real(ZETA,kind=4), trim(cfilvor), &
            & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_vor), trim(cunit_vor), &
            & trim(cln_vor), 0.0)

       CALL P2D_T(idx_f6, idx_v6, nx, nz, nt, icount, real(vlon,kind=4), real(vht,kind=4),&
            & real(vtime,kind=4), real(DIV,kind=4), trim(cfildiv), &
            & trim(cvarlon), trim(cvarht), trim(cvartime), trim(cvarout_div), trim(cunit_div), &
            & trim(cln_div), 0.0)

    END IF
    !
    ! write diagnostic of the last time step
    IF ( ptcounter == rn_tmax ) THEN
       OPEN(unit=numdiag, File="diag_minmax.dat", STATUS="REPLACE", FORM="FORMATTED", ACCESS="SEQUENTIAL")

       WRITE(numdiag,*)
       WRITE(numdiag,*)'time to write the maximum and minimum of all fields'

       WRITE(numdiag,*)'----------------------------------------------------'
       WRITE(numdiag,*)'maximum value of each field :'
       WRITE(numdiag,*)
       WRITE(numdiag,*)"max_U           =", REAL(MAXVAL(U(:,:)), wp)
       WRITE(numdiag,*)"max_W           =", REAL(MAXVAL(W(:,:)), wp)
       WRITE(numdiag,*)"max_P_prime     =", REAL(MAXVAL(P(:,:)-PH(:,:)), wp)
       WRITE(numdiag,*)"max_THETA_prime =", REAL(MAXVAL(THETA(:,:)-TP(:,:)), wp)
       WRITE(numdiag,*)"max_ZETA        =", REAL(MAXVAL(ZETA(:,:)), wp)
       WRITE(numdiag,*)"max_DIV         =", REAL(MAXVAL(DIV(:,:)), wp)
       WRITE(numdiag,*)

       WRITE(numdiag,*)'----------------------------------------------------'
       WRITE(numdiag,*)'minimum value of each field :'
       WRITE(numdiag,*)
       WRITE(numdiag,*)"min_U           =", REAL(MINVAL(U(:,:)), wp)
       WRITE(numdiag,*)"min_W           =", REAL(MINVAL(W(:,:)), wp)
       WRITE(numdiag,*)"min_P_prime     =", REAL(MINVAL(P(:,:)-PH(:,:)), wp)
       WRITE(numdiag,*)"min_THETA_prime =", REAL(MINVAL(THETA(:,:)-TP(:,:)), wp)
       WRITE(numdiag,*)"min_ZETA        =", REAL(MINVAL(ZETA(:,:)), wp)
       WRITE(numdiag,*)"min_DIV         =", REAL(MINVAL(DIV(:,:)), wp)
    END IF

  END SUBROUTINE cycl_mck
  !!==================================================================
END MODULE dcn_cycle
