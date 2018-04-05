MODULE phys_cst
  !!==================================================================
  !!                 *** MODULE phys_cst ***
  !!      DEFINITION of constant variable for dry convection model
  !!==================================================================

  !!------------------------------------------------------------------
  !! rgamma_com : compute ratio of specific heats of air rcp/rcv
  !! read_lapsrt: read the lpaserate of basic state from namelist
  !!------------------------------------------------------------------
  USE par_kind      ! kind parameter of variables
  USE in_out_manager ! unit for namelist file 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: comp_rgamma, read_lapsrt

  REAL(wp), PUBLIC, PARAMETER ::  rpi  = 3.141592653589793    !: pi
  REAL(wp), PUBLIC, PARAMETER ::  grav = 9.80665_wp      !:gravity (m/s2)
  REAL(wp), PUBLIC, PARAMETER ::  rcp  = 1004.0_wp       !: specific_heat_air_pressure
  REAL(wp), PUBLIC, PARAMETER ::  rcv  = 717.0_wp        !: specific_heat_air_volume 
  REAL(wp), PUBLIC, PARAMETER ::  rps  = 1.0e5_wp           !: assumed_surface_pressure  
  REAL(wp), PUBLIC, PARAMETER ::  rd   = 287.04_wp       !: gas_constant_dry_air
  REAL(wp), PUBLIC, PARAMETER ::  rthetas = 300.0_wp     !: potential_temperature_surface 
  REAL(wp), PUBLIC ::  rgamma                            !: ratio_specific_heats_air
  REAL(wp), PUBLIC ::  rn_lapsrate                       !: lapserate_initial_state
  REAL(wp), PUBLIC ::  rn_nu                             !: numerical diffusivity constant
  LOGICAL,  PUBLIC ::  ln_neutral                        !: flag to check neutrality

CONTAINS

  SUBROUTINE comp_rgamma 
    !!--------------------------------------------------------------
    !!              *** ROUTINE rgamma ***   
    !! ** Purpose : compute ratio of specific heats of air
    !!--------------------------------------------------------------

    rgamma = rcp / rcv
    !
  END SUBROUTINE comp_rgamma

  SUBROUTINE read_lapsrt 
    !!--------------------------------------------------------------
    !!              *** ROUTINE read_lapsrt ***   
    !! ** Purpose : read the lpaserate of basic state from name list
    !!--------------------------------------------------------------
    NAMELIST/namstab/rn_lapsrate, ln_neutral, rn_nu
    !
    !!OPEN(unit = numnam, FILE='namelist', status="unknown")
    !
    REWIND(numnam)
    READ(numnam, namstab)
    !
    !!CLOSE(numnam)
    !
  END SUBROUTINE read_lapsrt

  !!==================================================================
END MODULE phys_cst
