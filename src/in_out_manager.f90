MODULE in_out_manager
  !!====================================================================================
  !!                       ***  MODULE  in_out_manager  ***
  !!               Declaration of logical units for all the files 
  !!                used as either in or out in dcn model
  !!====================================================================================
  USE par_kind
  
  IMPLICIT NONE
  PRIVATE

  !!-------------------------------------------------------------------------------------
  !!                                integer units
  !!-------------------------------------------------------------------------------------

  INTEGER, PUBLIC :: numscr  = 6       !: logical unit for screen
  INTEGER, PUBLIC :: numstp  = 10      !: logical unit for time.step
  INTEGER, PUBLIC :: numout  = 11      !: logical unit for dcn.output
  INTEGER, PUBLIC :: numnam  = 12      !: logical unit for namelist
  INTEGER, PUBLIC :: numu    = 13      !: logical unit for u variable
  INTEGER, PUBLIC :: numw    = 14      !: logical unit for w variable
  INTEGER, PUBLIC :: nump    = 15      !: logical unit for pressure 
  INTEGER, PUBLIC :: numpt   = 16      !: logical unit for potential_temperature 
  INTEGER, PUBLIC :: numvor  = 17      !: logical unit for vorticity
  INTEGER, PUBLIC :: numdiv  = 18      !: logical unit for divergence
  INTEGER, PUBLIC :: numdiag = 19      !: logical unit for diag_minmax
  !!--------------------------------------------------------------------------------------
  !!                       character for output file names  
  !!--------------------------------------------------------------------------------------
  CHARACTER(len=300), PUBLIC :: cn_drout
  CHARACTER(len=300), PUBLIC :: cfilu
  CHARACTER(len=300), PUBLIC :: cfilw
  CHARACTER(len=300), PUBLIC :: cfilpt
  CHARACTER(len=300), PUBLIC :: cfilp
  CHARACTER(len=300), PUBLIC :: cfilvor
  CHARACTER(len=300), PUBLIC :: cfildiv
  CHARACTER(len=30),  PUBLIC :: cvarht   = 'ht'
  CHARACTER(len=30),  PUBLIC :: cvarlon  = 'lon'
  CHARACTER(len=30),  PUBLIC :: cvartime = 'time'
  CHARACTER(len=30),  PUBLIC :: cunit_u, cunit_w, cunit_pt, cunit_p, &
       &cunit_vor, cunit_div
  CHARACTER(len=50),  PUBLIC :: cvarout_u, cvarout_w, cvarout_pt, cvarout_p, &
       &cvarout_vor, cvarout_div
  CHARACTER(len=300), PUBLIC :: cln_u, cln_w, cln_pt, cln_p, &
       & cln_vor, cln_div
  !!---------------------------------------------------------------------------------------
  !!                                integer units for netcdf
  !!======================================================================================
 
  INTEGER, PUBLIC :: nt   !: The number of time step to be written on netcdf file
  INTEGER, PUBLIC :: idx_f1, idx_f2, idx_f3, idx_f4, &
       idx_f5, idx_f6
  INTEGER, PUBLIC :: idx_v1, idx_v2, idx_v3, idx_v4, &
       idx_v5, idx_v6
  REAL(wp), PUBLIC, DIMENSION(:), ALLOCATABLE :: vtime
END MODULE in_out_manager
