MODULE dcn_opt
  !!======================================================================
  !!                    ***  MODULE dcn_opt    ***
  !! differential operator :  defining diferential operator used by dcn
  !!======================================================================

  USE par_kind
  IMPLICIT NONE
  PRIVATE

  INTERFACE optfx
     MODULE PROCEDURE forwardx
  END INTERFACE

  INTERFACE optbx
     MODULE PROCEDURE backwardx
  END INTERFACE

  INTERFACE optfz
     MODULE PROCEDURE forwardz
  END INTERFACE

  INTERFACE optbz
     MODULE PROCEDURE backwardz
  END INTERFACE

  INTERFACE optf
     MODULE PROCEDURE forward
  END INTERFACE

  INTERFACE optb
     MODULE PROCEDURE backward
  END INTERFACE

  INTERFACE divdisvor
     MODULE PROCEDURE divergence, dissipation 
  END INTERFACE
  
  INTERFACE vor
     MODULE PROCEDURE  vorticity
  END INTERFACE
  
  PUBLIC :: forwardx, backwardx
  PUBLIC :: forwardz, backwardz
  PUBLIC :: forward, backward
  PUBLIC :: divergence, dissipation, vorticity

CONTAINS

  SUBROUTINE forwardx ( pdx, pf )
    !!------------------------------------------------------------------------------
    !!                  ***  ROUTINE forwardx  ***
    !!
    !! ** Purpose : Calculate first derivative of function pf using implicit &
    !!             &forward differencing
    !!
    !! ** Method  : uing fourth order implicit MacCormack type forward &
    !!              &differencing 
    !!
    !!
    !! ** Reference :
    !!
    !!----------------------------------------------------------------------
    REAL(wp), INTENT(in )   :: pdx
    REAL(wp), DIMENSION(:), INTENT(inout) :: pf
    !!
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: zdff
    REAL(wp)                             :: zalpha
    INTEGER :: ni, ji

    !!---------------------------------------------------------------------
    ni = SIZE(pf,1) 
    ALLOCATE( zdff(ni) )

    zalpha = 1.0_wp / 2.0_wp - 1.0_wp / (2.0_wp * SQRT(3.0_wp))

    zdff(1) = (-25.0_wp / 12.0_wp + 17.0_wp / (12.0_wp * SQRT (3.0_wp))) * pf(1) + &
         & (4.0_wp - 25.0_wp / (6.0_wp * SQRT (3.0_wp)) ) * pf(2) - &
         &( 3.0_wp - 3.0_wp * SQRT (3.0_wp) / 2.0_wp ) * pf(3) + ( 4.0_wp / 3.0_wp - 13.0_wp/ &
         &( 6.0_wp * SQRT(3.0_wp)) ) * pf(4) - ( 1.0_wp / 4.0_wp - &
         & 5.0_wp / ( 12.0_wp * SQRT(3.0_wp) ) ) * pf(5)

    zdff(1) = zdff(1) / pdx

    zdff(ni) = ( 25.0_wp / 12.0_wp + 17.0_wp / (12.0_wp * SQRT(3.0_wp) ) ) * pf(ni)- &
         &( 4.0_wp + 25.0_wp / ( 6.0_wp * SQRT(3.0_wp) ) ) * pf(ni-1) + &
         &( 3.0_wp + 3.0_wp * SQRT(3.0_wp) / 2.0_wp ) * pf(ni-2) - ( 4.0_wp / 3.0_wp &
         & + 13.0_wp / ( 6.0_wp * SQRT(3.0_wp) ) ) * pf(ni-3) + &
         &( 1.0_wp / 4.0_wp + 5.0_wp / (12.0_wp * SQRT(3.0_wp) ) ) * pf(ni-4)

    zdff(ni) = zdff(ni) / pdx

    zdff(ni) = 0.0_wp

    DO ji = ni-1, 2, -1
       zdff(ji) = ( 1.0_wp / pdx * ( pf(ji+1) - pf(ji)) - zalpha &
            & * zdff(ji+1)) / ( 1.0_wp - zalpha )
    END DO

    pf(:) = zdff(:)

    DEALLOCATE ( zdff )
    !
  END SUBROUTINE forwardx


  SUBROUTINE backwardx ( pdx, pf )
    !!------------------------------------------------------------------------------
    !!                  ***  ROUTINE forwardx  ***
    !!
    !! ** Purpose : Calculate first derivative of function pf using implicit &
    !!             &forward differencing
    !!
    !! ** Method  : uing fourth order implicit MacCormack type forward &
    !!              &differencing 
    !!
    !!
    !! ** Reference :
    !!
    !!----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: pdx
    REAL(wp), DIMENSION(:), INTENT(inout) :: pf
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: zdbf
    REAL(wp)                             :: zalpha
    INTEGER :: ni, ji
    !!---------------------------------------------------------------------------
    ni = SIZE(pf,1)
    ALLOCATE( zdbf(ni) )

    zalpha = 1.0_wp / 2.0_wp - 1.0_wp / (2.0_wp * SQRT(3.0_wp))

    zdbf(1) = - ( 25.0_wp / 12.0_wp + 17.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(1) &
         & + ( 4.0_wp + 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(2)  - &
         & ( 3.0_wp + 3.0_wp * SQRT(3.0_wp) / 2.0_wp ) * pf(3) + ( 4.0_wp / 3.0_wp + &
         &13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(4) - ( 1.0_wp / 4.0_wp + &
         & 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(5)

    zdbf(1) = zdbf(1) / pdx
    zdbf(1) = 0.0_wp

    zdbf(ni) = ( 25.0_wp / 12.0_wp - 17.0_wp / (12.0_wp * SQRT(3.0_wp))) * pf(ni) - & 
         & ( 4.0_wp - 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-1) + &
         & ( 3.0_wp - 3.0_wp * SQRT ( 3.0_wp) / 2.0_wp ) * pf(ni-2) - &
         & (4.0_wp / 3.0_wp - 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-3) &
         & + ( 1.0_wp / 4.0_wp - 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(ni-4)

    zdbf(ni) = zdbf(ni) / pdx

    DO ji = 2, ni-1
       zdbf(ji) = ( 1.0_wp / pdx * ( pf(ji) - pf(ji-1)) - zalpha * zdbf(ji-1)) /&
            & ( 1.0_wp - zalpha )
    END DO

    pf(:) = zdbf(:)

    DEALLOCATE (zdbf)
    !
  END SUBROUTINE backwardx


  SUBROUTINE backwardz ( pdx, pf )
    !!------------------------------------------------------------------------------                                                                                                         

    REAL(wp), INTENT(in )   :: pdx
    REAL(wp), DIMENSION(:), INTENT(inout) :: pf
    !!                                                                                                                                                                                       
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: zdbf
    REAL(wp)                             :: zalpha
    INTEGER :: ni, ji
    !!----------------------------------------------------------------------------

    ni = SIZE(pf,1)
    ALLOCATE( zdbf(ni) )


    zalpha = 1.0_wp / 2.0_wp - 1.0_wp / (2.0_wp * SQRT(3.0_wp))

    zdbf(1) = - (25.0_wp / 12.0_wp + 17.0_wp / (12.0_wp * SQRT(3.0_wp))) * pf(1) + &
         & ( 4.0_wp + 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(2)-&
         & ( 3.0_wp + 3.0_wp * SQRT ( 3.0_wp) / 2.0_wp ) * pf(3) + &
         & ( 4.0_wp / 3.0_wp + 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(4) &
         &-( 1.0_wp / 4.0_wp + 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(5)

    zdbf(1) = zdbf(1)/pdx

    zdbf(ni) = ( 25.0_wp / 12.0_wp - 17.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(ni)- &
         &( 4.0_wp - 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-1) + &
         & ( 3.0_wp - 3.0_wp * SQRT( 3.0_wp) / 2.0_wp ) * pf(ni-2) - &
         &( 4.0_wp / 3.0_wp - 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-3) &
         & + ( 1.0_wp / 4.0_wp - 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(ni-4)

    zdbf(ni) = zdbf(ni) / pdx

    DO ji = 2, ni-1
       zdbf(ji) = ( 1.0_wp / pdx * ( pf(ji) - pf(ji-1)) - zalpha * &
            & zdbf(ji-1)) / ( 1.0_wp - zalpha )
    END DO

    pf(:) = zdbf(:)

    DEALLOCATE( zdbf)
    !
  END SUBROUTINE backwardz


  SUBROUTINE forwardz ( pdx, pf )
    !!------------------------------------------------------------------------------                                                                                                         
    !!------------------------------------------------------------------------------
    !!                  ***  ROUTINE forwardx  ***
    !!
    !! ** Purpose : Calculate first derivative of function pf using implicit &
    !!             &forward differencing
    !!
    !! ** Method  : uing fourth order implicit MacCormack type forward &
    !!              &differencing 
    !!
    !!
    !! ** Reference :
    !!
    !!----------------------------------------------------------------------

    !!------------------------------------------------------------------------------
    REAL(wp), INTENT(in )   :: pdx
    REAL(wp), DIMENSION(:), INTENT(inout) :: pf
    !!                                                                                                                                                                                       
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: zdff
    REAL(wp)                             :: zalpha
    INTEGER :: ni, ji
    !!-----------------------------------------------------------------------------
    ni = SIZE(pf,1)
    ALLOCATE( zdff(ni) )

    zalpha = 1.0_wp / 2.0_wp - 1.0_wp / ( 2.0_wp * SQRT(3.0_wp) )

    zdff(1) = ( -25.0_wp / 12.0_wp + 17.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(1) + &
         & ( 4.0_wp - 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(2) - &
         & ( 3.0_wp - 3.0_wp * SQRT(3.0_wp) / 2.0_wp) * pf(3) + ( 4.0_wp / 3.0_wp - &
         & 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(4) - ( 1.0_wp / 4.0_wp - &
         & 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(5)

    zdff(1) = zdff(1) / pdx

    zdff(ni) = ( 25.0_wp / 12.0_wp + 17.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(ni) &
         & - ( 4.0_wp + 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-1) + &
         & ( 3.0_wp + 3.0_wp * SQRT (3.0_wp) / 2.0_wp) * pf(ni-2) - ( 4.0_wp / &
         & 3.0_wp + 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-3) + ( 1.0_wp / &
         4.0_wp + 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(ni-4)

    zdff(ni) = zdff(ni) / pdx

    DO ji = ni-1, 2, -1
       zdff(ji) = ( 1.0_wp / pdx * ( pf(ji+1) - pf(ji)) &
            & - zalpha * zdff(ji+1)) / ( 1.0_wp - zalpha )
    END DO

    pf(: ) = zdff(:)



  END SUBROUTINE forwardz




  SUBROUTINE forward (pdx, pf)
    !!------------------------------------------------------------------------------                                                                                                         
    !!------------------------------------------------------------------------------
    !!                  ***  ROUTINE forwardx  ***
    !!
    !! ** Purpose : Calculate first derivative of function pf using implicit &
    !!             &forward differencing
    !!
    !! ** Method  : uing fourth order implicit MacCormack type forward &
    !!              &differencing 
    !!
    !!
    !! ** Reference :
    !!
    !!----------------------------------------------------------------------

    !!------------------------------------------------------------------------------
    REAL(wp), INTENT(in )   :: pdx
    REAL(wp), DIMENSION(:), INTENT(inout) :: pf
    !!                                                                                                                                                                                       
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: zdff
    REAL(wp)                             :: zalpha
    INTEGER :: ni, ji
    !!-----------------------------------------------------------------------------
    ni = SIZE(pf,1)
    ALLOCATE( zdff(ni) )

    zalpha = 1.0_wp / 2.0_wp-1.0_wp/(2.0_wp*SQRT(3.0_wp))

    zdff(1) =(-25.0_wp/12.0_wp+17.0_wp/(12.0_wp*SQRT(3.0_wp)))*pf(1)+(4.0_wp-25.0_wp/(6.0_wp*SQRT(3.0_wp)))*pf(2)-&
         &(3.0_wp-3.0_wp*SQRT(3.0_wp)/2.0_wp)*pf(3)+(4.0_wp/3.0_wp-13.0_wp/(6.0_wp*SQRT(3.0_wp)))*pf(4)-(1.0_wp/4.0_wp-&
         &5.0_wp/(12.0_wp*SQRT(3.0_wp)))*pf(5)

    zdff(1)=zdff(1)/pdx

    zdff(ni) = ( 25.0_wp / 12.0_wp + 17.0_wp / ( 12.0_wp * SQRT ( 3.0_wp ) ) ) * pf(ni) &
         &- (4.0_wp + 25.0_wp / ( 6.0_wp * SQRT (3.0_wp))) *pf(ni-1) + &
         &( 3.0_wp + 3.0_wp * SQRT(3.0_wp) / 2.0_wp ) * pf(ni-2) - ( 4.0_wp / 3.0_wp &
         & + 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-3) + (1.0_wp / 4.0_wp + &
         & 5.0_wp /( 12.0_wp * SQRT(3.0_wp))) * pf(ni-4)

    zdff(ni) = zdff(ni) / pdx

    zdff(ni) = 0.0_wp

    DO ji = ni-1, 2, -1
       zdff(ji) = (1.0_wp / pdx * ( pf(ji+1) - pf(ji)) - &
            & zalpha * zdff (ji+1)) / ( 1.0_wp -zalpha )
    END DO

    pf(:)=zdff(:)



  END SUBROUTINE forward


  SUBROUTINE backward (pdx, pf)
    !!------------------------------------------------------------------------------                                                                                                         
   !!                  ***  ROUTINE backwardz  ***                                                                                                                                           
    
    !! ** Purpose : Calculate first order derivative of function pf using implicit &                                                                                                         
    !!             &forward differencing                                                                                                                                                    

    !! ** Method  : uing fourth order implicit MacCormack type forward &                                                                                                                     


    REAL(wp), INTENT(in )   :: pdx
    REAL(wp), DIMENSION(:), INTENT(inout) :: pf
    !!                                                                                                                                                                                       
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: zdbf
    REAL(wp)                             :: zalpha
    INTEGER :: ni, ji
    !!-----------------------------------------------------------------------------
    ni = SIZE(pf,1)
    ALLOCATE( zdbf(ni) )


    zalpha = 1.0_wp / 2.0_wp - 1.0_wp / ( 2.0_wp * SQRT (3.0_wp) )

    zdbf(1) = - ( 25.0_wp / 12.0_wp + 17.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(1)&
         & + ( 4.0_wp + 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(2) - &
         & ( 3.0_wp + 3.0_wp * SQRT(3.0_wp) / 2.0_wp ) * pf(3) + (4.0_wp / 3.0_wp  + 13.0_wp / &
         & ( 6.0_wp * SQRT(3.0_wp))) * pf(4) - (1.0_wp / 4.0_wp  + &
         & 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(5)

    zdbf(1) = zdbf(1) / pdx

    zdbf(1) = 0.0_wp

    zdbf(ni) = ( 25.0_wp / 12.0_wp - 17.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * &
         & pf(ni) - ( 4.0_wp - 25.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-1) + &
         & ( 3.0_wp - 3.0_wp * SQRT(3.0_wp) / 2.0_wp ) * pf(ni-2) - ( 4.0_wp / 3.0_wp - &
         & 13.0_wp / ( 6.0_wp * SQRT(3.0_wp))) * pf(ni-3) + ( 1.0_wp / 4.0_wp - &
         & 5.0_wp / ( 12.0_wp * SQRT(3.0_wp))) * pf(ni-4)

    zdbf(ni) = zdbf(ni) / pdx

    DO ji = 2, ni-1
       zdbf(ji) = ( 1.0_wp / pdx * ( pf(ji) - pf(ji-1)) - &
            & zalpha * zdbf(ji-1)) / ( 1.0_wp - zalpha )
    END DO

    pf(:) = zdbf(:)



  END SUBROUTINE backward





  SUBROUTINE divergence (pdx, pdz, pu, pw, pdiv)

    REAL(wp), DIMENSION(:,:), INTENT(in)  :: pu, pw
    REAL(wp), DIMENSION(:,:), INTENT(out) :: pdiv
    REAL(wp), INTENT(in) :: pdx, pdz
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: zx, zz
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: uxb, uxf
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: wzf, wzb 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ux, wz   
    INTEGER :: ni, nj, ji, jj

    ni = SIZE( pu, 1 )
    nj = SIZE( pu, 2 )

    ALLOCATE ( zx(ni), zz(nj) )
    ALLOCATE( uxb(ni,nj), uxf(ni,nj))
    ALLOCATE( wzf(ni,nj), wzb(ni,nj))
    ALLOCATE( ux(ni,nj), wz(ni,nj))

    DO jj = 1, nj

       zx(:) = pu(:,jj)

       CALL forward ( pdx, zx )

       uxf(:,jj) = zx(:)

    END DO
    DO jj = 1, nj

       zx(:) = pu(:,jj)

       CALL backward ( pdx, zx )

       uxb(:,jj) = zx(:)

    END DO


    ux(:,:) = ( uxf(:,:) + uxb(:,:) ) / 2.0_wp


    DO ji = 1, ni

       zz(:) = pw(ji,:)

       CALL forward(  pdz, zz )

       wzf(ji,:) = zz(:)

    END DO
    DO ji = 1, ni

       zz(:) = pw(ji,:)

       CALL backward ( pdz, zz )

       wzb(ji,:) = zz(:)

    END DO

    wz(:,:) = ( wzf(:,:) + wzb(:,:) ) / 2.0_wp



    pdiv(:,:) = ux(:,:) + wz(:,:)

    DEALLOCATE ( zx, zz )
    DEALLOCATE( uxb, uxf)
    DEALLOCATE( wzf, wzb)
    DEALLOCATE( ux, wz)



  END SUBROUTINE divergence


  SUBROUTINE dissipation ( pdx, pdz, pu, prho, pnu, pdis )

    REAL(wp), DIMENSION(:,:), INTENT(in)  :: pu, prho
    REAL(wp), DIMENSION(:,:), INTENT(out) :: pdis
    REAL(wp), INTENT(in) ::   pdx, pdz, pnu
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: zx, zz
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: uxb, uxf
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: uzf, uzb 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ux, uz   
    INTEGER :: ni, nj, ji, jj

    ni = SIZE( pu, 1 )
    nj = SIZE( pu, 2 )

    ALLOCATE ( zx(ni), zz(nj) )
    ALLOCATE( uxb(ni,nj), uxf(ni,nj))
    ALLOCATE( uzf(ni,nj), uzb(ni,nj))
    ALLOCATE( ux(ni,nj), uz(ni,nj))



    DO jj = 1, nj

       zx(:) = pu(:,jj)

       CALL forward ( pdx, zx )

       uxf(:,jj) = zx(:)

    END DO

    DO jj = 1, nj

       zx(:) = pu(:,jj)

       CALL backward ( pdx, zx)
       DO ji = 1, ni
          uxb(:,jj) = zx(:)
       END DO
    END DO


    ux(:,:) = ( uxf(:,:) + uxb(:,:) ) / 2.0_wp


    DO ji = 1, ni

       zz(:) = pu(ji,:)

       CALL forward ( pdz, zz )
       DO jj = 1, nj
          uzf(ji,:) = zz(:)
       END DO
    END DO
    DO ji = 1, ni

       zz(:) = pu(ji,:)

       CALL backward ( pdz, zz )
       DO jj = 1, nj
          uzb(ji,:) = zz(:)
       END DO
    END DO

    uz(:,:) = ( uzf(:,:) + uzb(:,:) ) / 2.0_wp



    ux(:,:) = ux(:,:) * prho(:,:)
    uz(:,:) = uz(:,:) * prho(:,:)

    CALL divergence(  pdx, pdz, ux, uz, pdis)

    pdis(:,:) = pnu * pdis(:,:)

    DEALLOCATE ( zx, zz )
    DEALLOCATE( uxb, uxf)
    DEALLOCATE( uzf, uzb)
    DEALLOCATE( ux, uz)

  END SUBROUTINE dissipation

  SUBROUTINE vorticity ( pdx, pdz, pu, pw, pvory)

    REAL(wp), DIMENSION(:,:), INTENT(in)  :: pu, pw
    REAL(wp), DIMENSION(:,:), INTENT(out) :: pvory
    REAL(wp), INTENT(in) :: pdx, pdz
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: zx, zz
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: uzb, uzf
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: wxf, wxb 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: uz, wx   
    INTEGER :: ni, nj, ji, jj
    ni = SIZE( pu, 1 )
    nj = SIZE( pu, 2 )

    ALLOCATE ( zx(ni), zz(nj) )
    ALLOCATE( wxb(ni,nj), wxf(ni,nj))
    ALLOCATE( uzf(ni,nj), uzb(ni,nj))
    ALLOCATE( uz(ni,nj), wx(ni,nj))

    DO jj = 1, nj

       zx(:)= pw(:,jj)

       CALL forward ( pdx, zx ) 

       wxf(:,jj) = zx(:)

    END DO

    DO jj = 1, nj

       zx(:) = pw(:,jj)

       CALL backward ( pdx, zx )

       wxb(:,jj) = zx(:)

    END DO


    wx(:,:) = ( wxf(:,:) + wxb(:,:) ) / 2.0_wp


    DO ji = 1, ni

       zz(:) = pu(ji,:)

       CALL forward (pdz, zz)

       uzf(ji,:) = zz(:)

    END DO

    DO ji = 1, ni

       zz(:) = pu(ji,:)

       CALL backward ( pdz, zz )

       uzb(ji,:) = zz(:)

    END DO


    uz(:,:) = ( uzf(:,:) + uzb(:,:) )/ 2.0_wp



    pvory(:,:) = uz(:,:) - wx(:,:)

  END SUBROUTINE vorticity

END MODULE dcn_opt
