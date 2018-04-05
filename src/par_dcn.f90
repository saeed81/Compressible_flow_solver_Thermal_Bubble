MODULE par_dcn
  !!======================================================================
  !!                        ***  par_oce  ***
  !! Atmosphere :   set the atmosphere parameters
  !!
  !! References :
  !! 1-Mendenz-Nunez, L. R., 1990: Application of the MacCormack Scheme to Atmospheric Nonhydrostatic Models.&
  !!   Mon. Wea. Rev., vol. 122, 984-1000.
  !!
  !! 2-Droegemeier, K. K., and R. B. Wilhemson, 1987: Numerical simulation of thunderstorm outflow dynamics. Part I: & 
  !!  &Outflow sensitivity experiments and turbulent dynamics. J. Atmos. Sci., vol. 44, 1180-1210.
  !!
  !! 3-Straka JM, Wilhelmson RB, Wicker LJ, Anderson JR, Droegemeier KK. Numerical Solutions of A &
  !! & Non-linear Density Current: A Benchmark Solution and Comparisons. International Journal for &
  !!   Numerical Methods in Fluids 1993; 17:1–22.
  !!
  !!======================================================================
  USE par_kind      !kind parameter
  USE in_out_manager
  
  IMPLICIT NONE
  PUBLIC :: get_cexper, get_config

  !! Atmosphere Domain sizes
  !! number of grid points
  !! ------------------

  INTEGER, PUBLIC, PARAMETER  ::      &  !:
       & nx = 501,                   & !: number of grid points in x direction 
       & nz = 201                       !: number of grid points in z direction

  !! Atmosphere Domain sizes
  !! length and width of domain
  !!---------------------

  !! 4 widely used test cases are as follows ::
  !! 1-straka                        !:straka lx = 25600.0_wp, lz = 6400.0_wp
  !! 2-density current               !:dnctr  lx = 25000.0_wp, lz = 10000.0_wp
  !! 3-warm bubblw with bottom rigid !:wbrbb  lx = 40000.0_wp, lz = 15000.0_wp
  !! 4-warm bubble with 4 rigid      !:wbr4b  lx = 3200.0_wp, lz = 4000.0_wp
  !! --------------------

  REAL(wp), PUBLIC ::       &
       & rn_lx ,             &  !: length of domain in x direction
       & rn_lz                   !: length of domain in z direction

  !!Atmospheric domain size
  !!grid spacing
  !!----------------------

  REAL(wp), PUBLIC ::       &
       & dx ,             &  !: grid spacing of domain in x direction
       & dz                  !: grid spacing of domain in z direction

  !!------------------------
  !!The name of numerical experiment
  !!wbrbb  !: warm bubble in a rectangular domain with rigid bottom boundary                       
  !!wbr4b  !: warm bubble in a rectangular domain with 4 rigid boundary
  !!straka !: straka test case
  !!dnctr  !: density current
  !!------------------------ 
  CHARACTER(len=10), PUBLIC :: &
       & cn_exper         !: the name of experience
  !!----------------------
  
  !!----------------------
  !!flag for rigid boundary 
  !!True means rigid otherwise open boundary
  !!straka  lrb=.TRUE., lru=, lrl=,lrr=
  !!dnctr   lrb=,lru=,lrl=,lrr=
  !!wbrbb   lrb=,lru=,lrl=,lrr=
  !!wbr4b   lrb=,lru=,lrl=,lrr=

  LOGICAL, PUBLIC :: &
       & ln_rb          !:flag for rigid bottom boundary

  LOGICAL, PUBLIC :: &   
       &ln_ru          !:flag for rigid upper boundary

  LOGICAL, PUBLIC :: &
       &ln_rl          !:flag for rigid left boundary 

  LOGICAL, PUBLIC :: &
       &ln_rr          !:flag for rigid right boundary

CONTAINS 

  SUBROUTINE get_cexper

    NAMELIST/namcexper/cn_exper

    REWIND(numnam)
    READ(numnam,namcexper)

  END SUBROUTINE get_cexper

  SUBROUTINE get_config

    NAMELIST/namconfig/rn_lx, rn_lz, ln_rb, ln_ru &
         &, ln_rl, ln_rr

    REWIND(numnam)
    READ(numnam,namconfig)

    dx = rn_lx / (nx - 1)
    dz = rn_lz / (nz - 1)

  END SUBROUTINE get_config

  !!==========================================================================
END MODULE par_dcn
