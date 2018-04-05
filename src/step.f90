MODULE step
  USE par_kind
  USE in_out_manager

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: par_step

  REAL(wp), PUBLIC :: rn_dt
  REAL(wp), PUBLIC :: rn_tmax
  REAL(wp), PUBLIC :: rn_twrite
  INTEGER,  PUBLIC :: nn_tmax 
  INTEGER,  PUBLIC :: nn_twrite 

CONTAINS

  SUBROUTINE par_step

    NAMELIST/namstep/rn_dt, rn_tmax, rn_twrite, nn_tmax, nn_twrite

    REWIND(numnam)
    READ(numnam,namstep)

  END SUBROUTINE par_step

END MODULE step
