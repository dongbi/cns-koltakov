// Parameters class  
// num_total_nodes_(x,y,z) must be divisible by 16*num_cpu_(x,y,z)

#ifndef __PARAMETERS__
#define __PARAMETERS__

#include <cmath>
#include <ctime>
#include <string>
#include <stdlib.h>
#include "array_2d.h"
#include "parameter_file_parser.h"

typedef enum {FREE_SLIP, NO_SLIP} BC_Type;
typedef enum {UPWIND, MUSCL, SHARP} Limiter_Type;

template<class T=double>
class PARAMETERS
{
 public:
 bool periodic_in_x, periodic_in_y, periodic_in_z, 
      resolve_interface_in_y, stretch_in_x, stretch_in_y, stretch_in_z, 
      potential_energy, scalar_advection, read_grid_from_file, aggregate_data,
      save_fluxes, save_instant_velocity, save_pressure, sediment_advection, 
      turbulence, moving_grid, open_top, variable_fixed_depth, coriolis;
 int  restart_timestep, max_timestep, 
      save_data_timestep_period, print_timestep_period,
      mg_sub_levels, max_mg_iters, mg_max_smoothing_iters, 
      mg_smoothing_sub_iters, 
      halo_size, num_local_nodes_x, num_local_nodes_y,num_local_nodes_z,
      num_cpu_x, num_cpu_y, num_cpu_z, num_total_nodes_x, num_total_nodes_y, 
      num_total_nodes_z, time_step, num_scalars;
  int i_min_w_h, i_max_w_h, j_min_w_h, j_max_w_h, k_min_w_h, k_max_w_h,
      i_min, i_max, j_min, j_max, k_min, k_max;
  T x_min, x_max, y_min, y_max, z_min, z_max, x_length, y_length, z_length,
    domain_skew_angle, y_stretching_ratio;
  T mg_smoothing_converg_thresh, mg_tol_absolute_resid, mg_tol_error_resid, 
    mg_tol_relative_resid, max_cfl, critical_cfl;
  T time, delta_time, molecular_viscosity, molecular_diffusivity, g, pi, 
    omega, amp_p_grad, freq_p_grad;
  std::string output_dir, grid_filename;
  int argc; 
  char** argv;

  clock_t sim_time;
  float elapsed_time, total_time;

  BC_Type west_bc, east_bc, suth_bc, nrth_bc, back_bc, frnt_bc;
  Limiter_Type universal_limiter;

  ARRAY_2D<T> *depth;
  ARRAY_2D<VECTOR_3D<T> > *lid_velocity, *bed_velocity;
  VECTOR_3D<T>  *pressure_gradient;

  PARAMETERS(int ac, char** av) : argc(ac),argv(av) 
  {
    parser = new PARAMETER_FILE_PARSER<T>("parameters.dat");
    parser->Parse_Parameter_File();
    Set_Parsable_Values();
    Set_Remaining_Parameters();
  }

  ~PARAMETERS()
   {delete parser;
    if(mg_sub_levels){
      delete [] num_subgrid_total_nodes_x; delete [] num_subgrid_total_nodes_y;
      delete [] num_subgrid_total_nodes_z;
      delete [] num_subgrid_local_nodes_x; delete [] num_subgrid_local_nodes_y;
      delete [] num_subgrid_local_nodes_z;}
    if(lid_velocity) delete lid_velocity;
    if(bed_velocity) delete bed_velocity;
    if(depth) delete depth;
    if(pressure_gradient) delete pressure_gradient;}

  void Set_Lid_Velocity(const VECTOR_3D<T>& v);
  void Set_Depth_Based_On_Node_Locations(const ARRAY_3D<VECTOR_3D<T> >& nodes);
 
 private:
  PARAMETER_FILE_PARSER<T> *parser;
  int *num_subgrid_total_nodes_x, *num_subgrid_total_nodes_y, 
      *num_subgrid_total_nodes_z, *num_subgrid_local_nodes_x, 
      *num_subgrid_local_nodes_y, *num_subgrid_local_nodes_z;  

  void Set_Parsable_Values();
  void Set_Remaining_Parameters();
  void Init_Depth_With_Random_Perturbation(const T pert_amplitude);
  void Init_Depth_With_Sinusoid_Perturbation_In_X_Direction(
	                    const T pert_amplitude, const int num_pert_periods, 
                            const ARRAY_3D<VECTOR_3D<T> >& nodes);
  void Init_Depth_With_Sloping_Bottom(
                            const ARRAY_3D<VECTOR_3D<T> >& nodes);

};
//*****************************************************************************
// Sets values of parameters based on 'parameters.dat'
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Set_Parsable_Values() {
  // domain size
  if(!parser->Get_Value("x_min",x_min)) x_min = 0.;
  if(!parser->Get_Value("x_max",x_max)) x_max = .2;
  if(!parser->Get_Value("y_min",y_min)) y_min = -.1; 
  if(!parser->Get_Value("y_max",y_max)) y_max = 0.;
  if(!parser->Get_Value("z_min",z_min)) z_min = 0.;
  if(!parser->Get_Value("z_max",z_max)) z_max = .1;
  // grid size
  if(!parser->Get_Value("num_total_nodes_x",num_total_nodes_x))
    num_total_nodes_x = 64;
  if(!parser->Get_Value("num_total_nodes_y",num_total_nodes_y))
    num_total_nodes_y = 32;
  if(!parser->Get_Value("num_total_nodes_z",num_total_nodes_z))
    num_total_nodes_z = 32;
  // number of CPUs
  if(!parser->Get_Value("num_cpu_x",num_cpu_x)) num_cpu_x = 4;
  if(!parser->Get_Value("num_cpu_y",num_cpu_y)) num_cpu_y = 2;
  if(!parser->Get_Value("num_cpu_z",num_cpu_z)) num_cpu_z = 2;
 
  if(!parser->Get_Value("delta_time",delta_time)) delta_time = .002;
  if(!parser->Get_Value("max_timestep",max_timestep)) max_timestep = 200000;
  if(!parser->Get_Value("save_data_timestep_period",save_data_timestep_period))
    save_data_timestep_period = 500;
  if(!parser->Get_Value("print_timestep_period",print_timestep_period))
    print_timestep_period = 1;
  if(!parser->Get_Value("molecular_viscosity",molecular_viscosity)) 
    molecular_viscosity = 1e-6;
  if(!parser->Get_Value("output_dir",output_dir)) output_dir = "./output/";
  if(!parser->Get_Value("restart_timestep",restart_timestep))restart_timestep=0;
}
//*****************************************************************************
// 1) Sets parameters not expected from 'parameters.dat'
// 2) Calculate auxiliary data structures
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Set_Remaining_Parameters(){
  // boolean parameters
  scalar_advection = true;
  num_scalars = 2; //1: rho only, 2: rho and passive scalar
  potential_energy = false; //true; //based on scalar
  sediment_advection = false;
  turbulence = false;
  moving_grid = false; //true;
  open_top = false;
  variable_fixed_depth = true; //sinusoidal bathymetry
  coriolis = false;  
  read_grid_from_file = false; //true;
  // boundary conditions
  periodic_in_x = false;
  periodic_in_y = false;
  periodic_in_z = true;
  resolve_interface_in_y = false;//true;// move nodes towards mid-depth
  stretch_in_x = false;   // move nodes towards the boundary
  stretch_in_y = false;
  stretch_in_z = false;
  y_stretching_ratio = (T)0.; //(T)1.03; // if =0, uniform in vertical
  west_bc = east_bc = NO_SLIP;
  //suth_bc = FREE_SLIP; // sloshing wave
  suth_bc = NO_SLIP; // lock-exchange
  nrth_bc = back_bc = frnt_bc = FREE_SLIP;
  universal_limiter = SHARP; //UPWIND; //MUSCL; //scalar convection

  // saving on disk
  //save_data_timestep_period = 500; //write on disk after each period
  save_fluxes = false; 
  save_instant_velocity = true; //save instantaneous velocity field
  save_pressure = false; //save pressure field
  aggregate_data = false; //save timeseries of any physical variables

  // multigrid
  mg_sub_levels = -1; //negative => optimal: calculated below
  max_mg_iters = 128;
  mg_max_smoothing_iters = 20; // total smoothing iters = iters*sub_iters
  mg_smoothing_sub_iters = 2; 
  mg_smoothing_converg_thresh = (T).6; //.006;
  mg_tol_absolute_resid = (T)1;
  mg_tol_error_resid = (T)1;
  mg_tol_relative_resid = (T)1e-5;

  time_step = 0;           // current time step
  time = (T)0;             // current time
  elapsed_time = 0;        // runtime per time step
  total_time = 0;          // total simulation runtime
  //delta_time = (T).002; //.01;(le.001)//.0001;(le)//.00001; // time increment
  //molecular_viscosity   = 1e-6; //(T).1;//1e-4;//5774e-12;//1e-6;//001;//.5;
  molecular_diffusivity = (T)0;//.5;
  pi = (T)3.14159265;
  g = (T)9.80665;
  omega = (T)0;      // Coriolis
  amp_p_grad = (T).5;
  freq_p_grad = (T).7852;
  max_cfl = (T)0;        // cfl accumulated for all time_steps
  critical_cfl = (T).98; // termination barrier for cfl
  lid_velocity = NULL;
  bed_velocity = NULL;
  depth = NULL;
  pressure_gradient = NULL;
  grid_filename = "grid_sloshing_wave.dat";//"grid_lock_exchange.dat";

  x_length = x_max - x_min;
  y_length = y_max - y_min; 
  z_length = z_max - z_min;
  domain_skew_angle = (T)90; //in degrees (if <90, domain is a parallelogram)

  // local nodes per proc (total nodes has to be divisible by num_cpu)
  assert(num_total_nodes_x % num_cpu_x == 0);
  assert(num_total_nodes_y % num_cpu_y == 0);
  assert(num_total_nodes_z % num_cpu_z == 0);
  num_local_nodes_x = num_total_nodes_x / num_cpu_x;
  num_local_nodes_y = num_total_nodes_y / num_cpu_y;
  num_local_nodes_z = num_total_nodes_z / num_cpu_z;

  // max allowable # of MG levels based on the number of local nodes in x/y/z
  if(mg_sub_levels<0){
    int mn = min( min(num_local_nodes_x,num_local_nodes_y), num_local_nodes_z );
    mg_sub_levels = log(mn)/log(2);
  }
  halo_size = 2; //number of halo cells per boundary for local grids
  i_min = j_min = k_min = 1;
  i_max=num_local_nodes_x; j_max=num_local_nodes_y; k_max=num_local_nodes_z;
  i_min_w_h = i_min-halo_size; i_max_w_h = i_max+halo_size,
  j_min_w_h = j_min-halo_size; j_max_w_h = j_max+halo_size;
  k_min_w_h = k_min-halo_size; k_max_w_h = k_max+halo_size;

  // set pressure gradient to drive the flow
  //pressure_gradient = new VECTOR_3D<T>(25e-5,0,0);

  //Set_Lid_Velocity(VECTOR_3D<T>(16,0,0));

  // setup structures for multigrid sublevels
  if(mg_sub_levels) {  
    //check if there are enough nodes for all sublevels (i.e., div 2^num_levs)
    int divisor = pow(2, mg_sub_levels);
    assert(num_local_nodes_x % divisor == 0);
    assert(num_local_nodes_y % divisor == 0);
    assert(num_local_nodes_z % divisor == 0);
    // setting up structures for multigrid subgrids
    num_subgrid_total_nodes_x = new int[mg_sub_levels];
    num_subgrid_total_nodes_y = new int[mg_sub_levels];
    num_subgrid_total_nodes_z = new int[mg_sub_levels];
    // subgrid nodes for the local cpu
    num_subgrid_local_nodes_x = new int[mg_sub_levels];
    num_subgrid_local_nodes_y = new int[mg_sub_levels];
    num_subgrid_local_nodes_z = new int[mg_sub_levels];

    for(int n = 0; n < mg_sub_levels; n++) {
      num_subgrid_total_nodes_x[n] = num_total_nodes_x / pow(2,n+1);
      num_subgrid_total_nodes_y[n] = num_total_nodes_y / pow(2,n+1);
      num_subgrid_total_nodes_z[n] = num_total_nodes_z / pow(2,n+1);  
      num_subgrid_local_nodes_x[n] = num_subgrid_total_nodes_x[n] / num_cpu_x;
      num_subgrid_local_nodes_y[n] = num_subgrid_total_nodes_y[n] / num_cpu_y;
      num_subgrid_local_nodes_z[n] = num_subgrid_total_nodes_z[n] / num_cpu_z;
    }
  }else
    num_subgrid_total_nodes_x = num_subgrid_total_nodes_y = 
    num_subgrid_total_nodes_z = num_subgrid_local_nodes_x =  
    num_subgrid_local_nodes_y = num_subgrid_local_nodes_z = NULL;
}
//*****************************************************************************
// Sets velocity of the top boundary: used for the lid-driven cavity flows
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Set_Lid_Velocity(const VECTOR_3D<T>& v)
{
  if(lid_velocity) delete lid_velocity;
  lid_velocity = new ARRAY_2D<VECTOR_3D<T> >(i_min,i_max,k_min,k_max,halo_size);
  for(int i=i_min_w_h; i<=i_max_w_h; i++)
    for(int k=k_min_w_h; k<=k_max_w_h; k++)
      (*lid_velocity)(i,k) = v;
}
//*****************************************************************************
// Callback function (called by the CURVILINEAR_GRID class)
// for setting depth based on location of grid nodes
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Set_Depth_Based_On_Node_Locations(
                                          const ARRAY_3D<VECTOR_3D<T> >& nodes)
{
  if(variable_fixed_depth)
   Init_Depth_With_Sloping_Bottom(nodes);
   //Init_Depth_With_Sinusoid_Perturbation_In_X_Direction(y_length/10.,2.,nodes);
}
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Init_Depth_With_Random_Perturbation(const T pert_amplitude)
{
  if(depth) delete depth;
  depth = new ARRAY_2D<T>(i_min, i_max, k_min, k_max, halo_size);
  T init_depth = y_length - pert_amplitude; //perturbation:[-pert_amp;+pert_amp]
  int i_min_w_h = 1-halo_size, i_max_w_h = num_local_nodes_x+halo_size,
      j_min_w_h = 1-halo_size, j_max_w_h = num_local_nodes_y+halo_size;
      k_min_w_h = 1-halo_size, k_max_w_h = num_local_nodes_z+halo_size;
  srand(time(NULL));

  for(int i=i_min_w_h; i<=i_max_w_h; i++)
    for(int k=k_min_w_h; k<=k_max_w_h; k++)
      (*depth)(i,k) = init_depth + 2*pert_amplitude * (rand() / (T)RAND_MAX);
}
//*****************************************************************************
// Sets depth(i,k) based on node locations for non-uniform grids, nodes:[0,1]
// Sinusoidal bottom
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Init_Depth_With_Sinusoid_Perturbation_In_X_Direction(
	                    const T pert_amplitude, const int num_pert_periods, 
	                    const ARRAY_3D<VECTOR_3D<T> >& nodes)
{
  if(depth) delete depth;
  depth = new ARRAY_2D<T>(i_min, i_max, k_min, k_max, halo_size);
  T period = (T)2 * pi * (T)num_pert_periods,
    shift = 0.;//(T).5/(T)num_total_nodes_x;
  //cout<<"nodes(i,j_min,k).x-shift="<<nodes(1,1,1).x-shift<<endl;
  for (int i=i_min_w_h; i<=i_max_w_h; i++)
    for (int k=k_min_w_h; k<=k_max_w_h; k++)
      (*depth)(i,k) = y_length - pert_amplitude +
	pert_amplitude * sin(period*(nodes(i,j_min,k).x-shift));
}
//*****************************************************************************
// Sets depth(i,k) based on node locations for non-uniform grids, nodes:[0,1]
// Sloping bottom
//*****************************************************************************
template<class T>
void PARAMETERS<T>::Init_Depth_With_Sloping_Bottom(
	                    const ARRAY_3D<VECTOR_3D<T> >& nodes)
{
  if(depth) delete depth;
  depth = new ARRAY_2D<T>(i_min, i_max, k_min, k_max, halo_size);

  T x_s = 1.; //x position at beginning of slope
  T node_s = x_s/x_length; //node at beginning of slope

  T rise = .5*y_length; //in physical space
  T run = x_length - x_s; //in physical space
  T slope = rise/(run/x_length); //in x node coordinates

  for (int i=i_min_w_h; i<=i_max_w_h; i++) {
    for (int k=k_min_w_h; k<=k_max_w_h; k++) {
      if (nodes(i,j_min,k).x < node_s)
        (*depth)(i,k) = y_length; 
      else
        (*depth)(i,k) = (y_length + slope*node_s) - slope*nodes(i,j_min,k).x;
      //(*depth)(i,k) = .25 - .25*tanh((nodes(i,j_min,k).x-1.));
      //(*depth)(i,k) = y_length - slope*nodes(i,j_min,k).x;
      //(*depth)(i,k) = y_length;
    }
  }
}
//*****************************************************************************
#endif
