MODULE dcn_vctr
  !!==================================================================
  !!                 *** MODULE dcn_vctr ***
  !!      explicit shape declaration  of all the arrays used for dcn
  !!==================================================================
  USE par_dcn

  IMPLICIT NONE
  PUBLIC
  
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: V  !:vector field  (rho, rho*u, rho*w, rho*theta)   
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: E  !:horizontal flux vector (u*rho, u*rho*u+p, u*rho*w, u*rho*theta)
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: F  !:vertical flux vector (w*rho, w*rho*u, w*rho*w, w*rho*theta)
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: H  !:buoyancy term  (0.0, 0.0, -rho*g, 0.0)
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: PP !:pressure term  (0.0, 0.0, p, 0.0)
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: Q  !:vector field varibale (rho, rho*u, rho*w, rho*theta)in
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: EIN !:horizontal flux vetor (u*rho, u*rho*u+p, u*rho*w, u*rho*theta)in
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: FIN !:vertical flux vector  (w*rho, w*rho*u, w*rho*w, w*rho*theta)in
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: HIN !:buoyancy term       (0.0, 0.0, -rho*g, 0.0)in
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: PPIN !:pressure term for w (0.0, 0.0, pin, 0.0)
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: DIS  !:dissipation term     (
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: DE   !:forward or backward differencing flux vector E, EIN    
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: DF   !:forward or backward differencing flux vector F, FIN   
  REAL(wp), PUBLIC, DIMENSION ( nx, nz, 4 ) :: DPP  !:forward or backward differencing for pressure term   
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: RHO  !:density                      [kg/m3]   
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: U    !:i-horizontal velocity        [m/s]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: W    !:vertical   velocity          [m/s]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: P    !:pressure                     [Pa]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: THETA !:Ptential temperature        [K]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: RHOTHTA !:rho*theta               [kg/m3 K]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: RHOIN !:density in                  [kg/m3]   
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: UIN   !:i-horizontal velocity in    [m/s]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: WIN   !:vertical velocity           [m/s]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: THETAIN !:Potential temperature in   [K]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: PIN     !:pressure in               [Pa]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: ZETA    !:vorticy in y direction    [s-1]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: DIV     !:divergence                [s-1]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: EXNER   !:exner function            [same as CP]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: D       !:used for dissipation operator
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: PH      !:initial hydrostatic pressure [Pa]
  REAL(wp), PUBLIC, DIMENSION ( nx, nz )    :: TP      !: initial potential temperature [K]
  REAL(wp), PUBLIC, DIMENSION ( nz )        :: vht     !: vector for height variable
  REAL(wp), PUBLIC, DIMENSION ( nx )        :: vlon    !: vector for horizontal x variable
  REAL(wp), PUBLIC, DIMENSION ( nx )        :: GX      !: used for differential opts    
  REAL(wp), PUBLIC, DIMENSION ( nz )        :: GZ      !: used for differential opts    
  REAL(wp), PUBLIC, DIMENSION ( nz )        :: RB      !: used for differential opts   
  REAL(wp), PUBLIC, DIMENSION ( nz )        :: RF      !: used for differential opts    
  REAL(wp), PUBLIC, DIMENSION ( nz )        :: ALPHA   !: used for hydrostatic term  [ non dimensional ]
  
  !!=============================================================================
END MODULE dcn_vctr

