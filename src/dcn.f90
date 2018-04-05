MODULE dcn
  !!===========================================================================================
  !!                       ***  MODULE dcn   ***
  !! Dry convection system : DCN dry convection system for 2 dimensional warm and cold bubbles
  !!===========================================================================================

  !!----------------------------------------------------------------------
  !!   dcn_model      : solve system of equations for dry convection
  !!   dcn_init       : initialization of the dcn model
  !!   dcn_closefile  : close remaining files
  !!----------------------------------------------------------------------
  USE dom_init   !initializing of dcn model (routine istate_init) 
  USE dcn_cycle  !time steeping part of dcn model (routine cycle_mck)
  USE iom        !input and output manager
  USE phys_cst
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: dcn_model ! called by model

CONTAINS

  SUBROUTINE dcn_model
    !!----------------------------------------------------------------------
    !!                     ***  ROUTINE dcn_model  ***
    !!
    !! ** Purpose :   dcn solves the 2-D nonhydrostatic compressible euler 
    !!                equation for dry air on the rectangular domain 
    !!
    !! ** Method  : - model general initialization
    !!              - launch the time-stepping (dcn_cycle routine)
    !!              - finalize the run by closing files 
    !!
    !! References : Falahat
    !!----------------------------------------------------------------------
    REAL(wp)        :: ztcounter ! dummy variable as a time counter
    INTEGER         :: istp        ! time step index
    !!------------------------------------------------------------------------
    !
    !                            !-----------------------!
    CALL dcn_init                !==  Initialisations  ==!
    !                            ! dom_init + write_init ! 
    !                            !-----------------------!
    !
    CALL FLUSH(numout)
    CALL FLUSH(numscr)
    !
    ztcounter = 0.0_wp
    istp = 0
    !
    !
    !                            !-----------------------!
    !                            !==   time stepping   ==!
    !                            !-----------------------!
    !add a flush of numout in case of problem in the init part...
    !!
    DO WHILE( ztcounter < rn_tmax )

       !! Cycle 1                                                                 
       !! ~~~~~~~                                                                 
       CALL cycl_mck ( istp, ztcounter, 1 )

       !! Cycle 2                                                                 
       !! ~~~~~~~                                                                 
       CALL cycl_mck ( istp, ztcounter, 2 )

       !! Cycle 3                                                                 
       !! ~~~~~~~                                                                 
       CALL cycl_mck ( istp, ztcounter, 3 )

       !! Cycle 4                                                                 
       !! ~~~~~~~                                                                 
       CALL cycl_mck ( istp, ztcounter, 4 )

    END DO
    !
    !                                 !--------------------------!
    !                                 !==  finalize the run    ==!
    !                                 !--------------------------!
    WRITE(numout,*) 'run stops at : ', istp * rn_dt 
    !
    CALL dcn_closefile 
    !  
  END SUBROUTINE dcn_model


  SUBROUTINE dcn_init
    !!----------------------------------------------------------------------
    !!                     ***  ROUTINE dcn_init  ***
    !!
    !! ** Purpose :   initialization of the dcn model
    !!
    !!----------------------------------------------------------------------
    !
    OPEN( unit = numout, FILE = 'dcn.output', FORM = "FORMATTED", STATUS='REPLACE') 
    !
    WRITE(numout,*)
    WRITE(numout,*) '        University Of Tehran             '
    WRITE(numout,*) '       Institute  Of Geophysics          '
    WRITE(numout,*) '         Dry Convection Model            '
    WRITE(numout,*) '          version 1.0  (2011)            '
    WRITE(numout,*) '              Authors :                  '
    WRITE(numout,*) '   Sarmad Ghader    sghader@ut.ac.ir     '
    WRITE(numout,*) '   Saeed  Falahat   saeed@misu.su.se     '
    WRITE(numout,*) '   Ali Bidokhti     bidokhti@ut.ac.ir    '
    WRITE(numout,*) '                                         '
    WRITE(numout,*) ' Initilization of dcn models starts      '
    !
    CALL istate_init
    !
    WRITE(numout,*) '  Initilization accomplished             '
    WRITE(numout,*) 'model is about to start time integration '
    WRITE(numout,*) '     the name of experiment              '
    WRITE(numout,*) '                                         '
    WRITE(numout,'(13x,a)') trim(cn_exper)
    WRITE(numout,*)
    WRITE(numout,*) 'rn_lapsrate', rn_lapsrate
    WRITE(numout,*) 'ln_neautral', ln_neutral
    WRITE(numout,*)
    WRITE(numout,*)'rn_nu = ', rn_nu
    WRITE(numout,*)'ln_rb = ', ln_rb
    WRITE(numout,*)'ln_ru = ', ln_ru
    WRITE(numout,*)'ln_rr = ', ln_rr
    WRITE(numout,*)'ln_rl = ', ln_rl
    WRITE(numout,*)
    WRITE(numout,*) '  the grid spacings  are as follows      '
    WRITE(numout,'(3x,"dx",1x,"=",F8.3,1x,"m")') dx 
    WRITE(numout,'(3x,"dz",1x,"=",F8.3,1x,"m")') dz 
    WRITE(numout,*)
    WRITE(numout,*) '  time parameters were set               '
    WRITE(numout,*) ' the time step of model                  '
    WRITE(numout,'(3x,"dt",1x,"=",F8.5,1x,"s")') rn_dt
    WRITE(numout,'(3x,"tmax",1x,"=",F9.3,1x,"s")') rn_tmax
    !
    OPEN(unit=numstp, FILE="time.step", STATUS="REPLACE", FORM="FORMATTED", ACCESS="SEQUENTIAL")
    !
  END SUBROUTINE dcn_init


  SUBROUTINE dcn_closefile
    !!-----------------------------------------------------------------------------
    !!                     ***  ROUTINE dcn_closefile  ***
    !!
    !! ** Purpose :   Close the files if they are open
    !!-----------------------------------------------------------------------------
    INTEGER :: ji          !dummy integer variable for loop
    INTEGER :: numq        !maximum number of quries for inquire command
    LOGICAL :: llopn, llok !logical variable for checking existenss and openness
    !!-----------------------------------------------------------------------------
    !
    numq = 100
    DO ji = 10, numq
       INQUIRE( unit = ji, exist = llok, opened = llopn)
       IF( llok .AND. llopn ) THEN
          CLOSE(ji)
       END IF
    END DO
    !      
  END SUBROUTINE dcn_closefile

  !!=========================================================================
END MODULE dcn
