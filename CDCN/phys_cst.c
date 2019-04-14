/*==================================================================
  !!                 *** MODULE phys_cst ***
  !!      DEFINITION of constant variable for dry convection model
  !!==================================================================

  !!------------------------------------------------------------------
  !! rgamma_com : compute ratio of specific heats of air rcp/rcv
  !! read_lapsrt: read the lpaserate of basic state from namelist
  !!------------------------------------------------------------------
*/
  
#define rpi     3.141592653589793  // pi
#define grav    9.806650           // gravity (m/s2)
#define rcp     1004.0             // specific_heat_air_pressure
#define rcv     717.0              // specific_heat_air_volume 
#define rps     1.0e50             // assumed_surface_pressure  
#define rd      287.040            // gas_constant_dry_air
#define rthetas 300.0              // potential_temperature_surface 
#define rgamma  (rcp / rcv)        // ratio_specific_heats_air
double rn_lapsrate = 0.0;          // lapserate_initial_state
double rn_nu       = 0.0;          // numerical diffusivity constant
int    ln_neutral  = 1;            // flag to check neutrality
