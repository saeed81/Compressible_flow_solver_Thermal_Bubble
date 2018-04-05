MODULE par_kind
  !!======================================================================
  !!                   ***  MODULE par_kind  ***
  !!         : define the kind of real for the whole model
  !!======================================================================
  !!
  IMPLICIT NONE
  PRIVATE


  ! Number model from which the SELECTED_*_KIND are requested:
  !             4 byte REAL       8 byte REAL
  ! CRAY:           -            precision = 13
  !                              exponent = 2465
  ! IEEE:      precision = 6     precision = 15
  !            exponent = 37     exponent = 307

  INTEGER, PUBLIC, PARAMETER ::        &    !: Floating point section
       sp = SELECTED_REAL_KIND( 6, 37),  &  !: single precision (real 4)
       dp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
       wp = dp                              !: working precision

  INTEGER, PUBLIC, PARAMETER ::        &    !: Integer section
       i4 = SELECTED_INT_KIND(9) ,       &  !: single precision (integer 4)
       i8 = SELECTED_INT_KIND(14)           !: double precision (integer 8)

  !!----------------------------------------------------------------------
END MODULE par_kind
