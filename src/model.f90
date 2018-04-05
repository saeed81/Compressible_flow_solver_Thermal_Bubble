PROGRAM model 
  !!----------------------------------------------------------------------     
  !!                    ***  PROGRAM model  ***                               
  !!   !! ** Purpose :   encapsulate the dcn model  
  !!----------------------------------------------------------------------      
  USE dcn    !dcn system   (dcn_model routine) 
  !!
  IMPLICIT NONE 
  !!
  CALL dcn_model   !dcn system
  !!
END PROGRAM model
