#include "navier_stokes_solver.h"

int main(int argc, char* argv[])
{
  NAVIER_STOKES_SOLVER<double> ns(argc,argv);
  // NS constructor parameters: 
  // (argc,argv,(physical domain boundaries), (# of grid points),(# of procs) )
  // lock-exchange example
  //NAVIER_STOKES_SOLVER<double> ns(argc,argv,0,.2,-.1,0,0,.1,128,64,64,4,2,2);
  //NAVIER_STOKES_SOLVER<double> ns(argc,argv,0,.8,0,.1,0,.1,256,64,8, 8,2,1); 
  //NAVIER_STOKES_SOLVER<double> ns(argc,argv,0,.8,0,.1,0,.1,512,128,8, 8,2,1);
  // sloshing wave example
  //NAVIER_STOKES_SOLVER<double> ns(argc,argv,0.,1.,-1.,0.,0.,1.,64,64,4,4,4,1);

  while( /*ns.No_NAN() &&*/ns.Check_CFL() && ns.Increment_Time_Step_Counter() ){  
    ns.Start_Simulation_Timer();
    ns.Move_Grid(); 
    ns.Predictor();
    ns.Enforce_Incompressibility();
    ns.Corrector();
    ns.Scalar_Solve();
    ns.Post_Process();
  }
  return 0;
}
