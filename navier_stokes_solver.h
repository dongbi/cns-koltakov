// Navier-Stokes solver main class
#ifndef __NAVIER_STOKES_SOLVER__
#define __NAVIER_STOKES_SOLVER__

#include "mpi_driver.h"
#include "curvilinear_moving_grid.h"
#include "convection.h"
#include "array_3d.h"
#include "vector_3d.h"
#include "turbulence.h"
#include "pressure.h"
#include "scalar.h"
#include "potential_energy.h"
#include "moving_grid_engine.h"
#include "data_aggregator.h"

template<class T=double>
class NAVIER_STOKES_SOLVER
{
 public:
  NAVIER_STOKES_SOLVER(int argc, char **argv);
 //NAVIER_STOKES_SOLVER(int argc, char **argv, T xmin, T xmax, T ymin, T ymax, 
 //           T zmin, T zmax, int num_nodes_x, int num_nodes_y, int num_nodes_z,
 //             int num_cpu_x, int num_cpu_y, int num_cpu_z);
  ~NAVIER_STOKES_SOLVER();

  void Predictor();
  void Enforce_Incompressibility() {pressure->Solve();}
  void Corrector();
  void Scalar_Solve();

  void Move_Grid() 
  { if(parameters->moving_grid) moving_grid_engine->Move_Grid_Nodes(); }

  T Background_Potential_Energy();
  T Calculate_CFL();
  T Total_Kinetic_Energy();

  bool Increment_Time_Step_Counter();
  bool Check_CFL();
  bool No_NAN();
  void Set_Initial_Conditions();
  void Set_Progressive_Wave_BC(T time);
  void Post_Process();
  void Start_Simulation_Timer();

 private:
  void Linear_Extrapolate_Into_Halo_Regions(ARRAY_3D<VECTOR_3D<T> >& u);
  void Enforce_Velocity_BC(ARRAY_3D<VECTOR_3D<T> >& u);

  void Add_Pressure_Gradient_Term(ARRAY_3D<VECTOR_3D<T> >& RHS);

  void Save_Simulation_Data();
  int Save_Binary_Simulation_Parameters();
  int Save_Binary_Simulation_Data();
  int Save_Binary_Simulation_Data(string a_name, ARRAY_3D<T>& a);
  int Save_Binary_Simulation_Data(string a_name, ARRAY_3D<VECTOR_3D<T> >& a);
  int Save_Simulation_Data_For_Restart();
  int Load_Simulation_Data_For_Restart(int restart_ts);

  PARAMETERS<T>* parameters;
  MPI_DRIVER<T>* mpi_driver;
  CURVILINEAR_GRID<T>* grid;
  CONVECTION<T>* convection1, *convection2;
  TURBULENCE<T>* turbulence;
  PRESSURE<T>* pressure;
  SCALAR<T>* scalar1, *scalar2;
  POTENTIAL_ENERGY<T>* potential_energy;
  DATA_AGGREGATOR<T>* data_aggregator;
  MOVING_GRID_ENGINE<T>* moving_grid_engine;

  ARRAY_3D<VECTOR_3D<T> > *u, *RHS_for_AB;
  ARRAY_1D<ARRAY_3D<T>* > *phi;
  ARRAY_3D<T> *P; //, *rho;
  ARRAY_3D<T> *U_xi, *U_et, *U_zt; //velocities on faces
};
//*****************************************************************************
template<class T>
void NAVIER_STOKES_SOLVER<T>::Scalar_Solve(){
  if(parameters->scalar_advection){
    if(parameters->num_scalars == 1){
      scalar1->Update_RHS();
      scalar1->Solve();
    }
    else if(parameters->num_scalars == 2){
      scalar1->Update_RHS();
      scalar2->Update_RHS();
      scalar1->Solve();
      scalar2->Solve();
    }
    if(parameters->potential_energy){
      T E_background = potential_energy->Calculate(); // executed by all
      T E_kinetic = Total_Kinetic_Energy();
      if(!mpi_driver->my_rank){
        cout<<"E_background = "<<E_background<<endl;
        cout<<"E_kinetic = "<<E_kinetic<<endl;
      }
      potential_energy->Write_To_Disk();
    }
  }
}
//*****************************************************************************
template<class T>
void NAVIER_STOKES_SOLVER<T>::Post_Process()
{
  // save data ever save_data_timestep_period
  if(parameters->time_step % parameters->save_data_timestep_period == 0){
    //Save_Simulation_Data();  
    Save_Binary_Simulation_Data();
    //Save_Simulation_Data_For_Restart();
  }

  if(parameters->aggregate_data) data_aggregator->Aggregate();
}
//*****************************************************************************
template<class T>
void NAVIER_STOKES_SOLVER<T>::Start_Simulation_Timer()
{
  if(parameters->time_step < parameters->max_timestep) 
    parameters->sim_time = clock();
}
//*****************************************************************************
template<class T>
bool NAVIER_STOKES_SOLVER<T>::Increment_Time_Step_Counter() 
{
  if(parameters->time_step < parameters->max_timestep) {
    parameters->time += parameters->delta_time;
    parameters->time_step++; 
    parameters->sim_time = clock()-parameters->sim_time;
    parameters->elapsed_time = ((float)parameters->sim_time)/CLOCKS_PER_SEC;
    parameters->total_time += parameters->elapsed_time;
    if(mpi_driver->my_rank == 0 && 
      parameters->time_step % parameters->print_timestep_period==0) {
      cout<< "Time step = " << parameters->time_step 
          << ", Time = "     << parameters->time << endl; 
      if(parameters->time_step > 1)
      cout<< "Previous time step took: " << parameters->elapsed_time << " s"
          << ", Total runtime: "   << parameters->total_time << " s" 
          << ", Time remaining ~ " << 
            (parameters->max_timestep-parameters->time_step)*parameters->elapsed_time/3600 
                                                                         << " hrs" << endl;
    }
    return true;
  }else{
    if(mpi_driver->my_rank == 0) cout<< "Reached MAX time step!" << endl;
    return false;
  }
}
//*****************************************************************************
template<class T>
bool NAVIER_STOKES_SOLVER<T>::Check_CFL() 
{
  T cfl = Calculate_CFL();
  if(mpi_driver->my_rank == 0 && 
       parameters->time_step % parameters->print_timestep_period==0) 
    cout<<"CFL = "<<cfl<<endl;
  if(cfl <= parameters->critical_cfl) {
    return true; 
  }else{
    if(mpi_driver->my_rank == 0) 
      cout<< "CFL=" << cfl << " is larger than allowed MAX="
          << parameters->critical_cfl << "!" << endl << "Aborting..." << endl;
    return false;
  }
}
//*****************************************************************************
template<class T>
bool NAVIER_STOKES_SOLVER<T>::No_NAN()
{
  T flag = 0.;
  if(u->Has_NAN_Values()){
    std::stringstream velocity_name;
    velocity_name << "velocity_nan_" << parameters->time_step;
    mpi_driver->Write_Global_Array_To_Disk(velocity_name.str(), *u);
    flag = 1.;
  }
  mpi_driver->Replace_With_Max_Value_Among_All_Procs(flag);
  if(flag != 0.) return false; else return true;
}
//*****************************************************************************
#endif
