//*****************************************************************************
// File:  curvilinear_moving_grid.h
// Class: CURVILINEAR_MOVING_GRID<T>
// Note:  descendant of the curvilinear_grid that has moving grid capabilities
//*****************************************************************************
#ifndef __CURVILINEAR_MOVING_GRID__
#define __CURVILINEAR_MOVING_GRID__

#include "curvilinear_grid.h"

template<class T=double>
class CURVILINEAR_MOVING_GRID : public CURVILINEAR_GRID<T>
{
 public:
  CURVILINEAR_MOVING_GRID(PARAMETERS<T> *params, MPI_DRIVER<T> *mpi);
  ~CURVILINEAR_MOVING_GRID();

  void Save_Grid_From_Previous_Timestep();
  void Update_Node_Velocities();
  void Update_Jacobian_Difference();

  void Write_Volume_Evolution_To_Disk();

  ARRAY_3D<T>  *Jacobian_diff,                     // cell volume change
               *U_grid_xi, *U_grid_et, *U_grid_zt, // grid node velocities 
               *prev_U_grid_xi, *prev_U_grid_et, *prev_U_grid_zt;
  ARRAY_3D<VECTOR_3D<T> > *prev_grid, *prev_XI, *prev_ET, *prev_ZT; 


  using CURVILINEAR_GRID<T>::grid;
  using CURVILINEAR_GRID<T>::XI_x; using CURVILINEAR_GRID<T>::XI_y;
  using CURVILINEAR_GRID<T>::XI_z;
  using CURVILINEAR_GRID<T>::ET_x; using CURVILINEAR_GRID<T>::ET_y;
  using CURVILINEAR_GRID<T>::ET_z;  
  using CURVILINEAR_GRID<T>::ZT_x; using CURVILINEAR_GRID<T>::ZT_y;
  using CURVILINEAR_GRID<T>::ZT_z;
  using CURVILINEAR_GRID<T>::inverse_Jacobian;
  using CURVILINEAR_GRID<T>::I_Min_With_Halo; using CURVILINEAR_GRID<T>::I_Min;
  using CURVILINEAR_GRID<T>::I_Max_With_Halo; using CURVILINEAR_GRID<T>::I_Max;
  using CURVILINEAR_GRID<T>::J_Min_With_Halo; using CURVILINEAR_GRID<T>::J_Min;
  using CURVILINEAR_GRID<T>::J_Max_With_Halo; using CURVILINEAR_GRID<T>::J_Max;
  using CURVILINEAR_GRID<T>::K_Min_With_Halo; using CURVILINEAR_GRID<T>::K_Min;
  using CURVILINEAR_GRID<T>::K_Max_With_Halo; using CURVILINEAR_GRID<T>::K_Max;
  using CURVILINEAR_GRID<T>::Halo_Size;

 private:
  void Update_Bed_And_Lid_Velocities(/*ARRAY_2D<T>& old_depth, */);
  void Set_Initial_Moving_Grid_Node_Positions();
  void Set_Node_Position_Example_1(T current_time);
  void Fix_New_Grid();

  using CURVILINEAR_GRID<T>::parameters;
  using CURVILINEAR_GRID<T>::mpi_driver;

  T* volume_evolution;
  int vol_array_size;
};
//*****************************************************************************
// Constructor
//*****************************************************************************
template<class T>
CURVILINEAR_MOVING_GRID<T>::CURVILINEAR_MOVING_GRID(PARAMETERS<T> *params, 
					            MPI_DRIVER<T> *mpi)
  : CURVILINEAR_GRID<T>(params->x_min, params->x_max,params->y_min, 
     params->y_max,params->z_min,params->z_max,params->i_min,params->i_max,
     params->j_min,params->j_max,params->k_min,params->k_max,params->halo_size,
     params->mg_sub_levels)
{
  mpi_driver = mpi; parameters = params;

  //local grid: [0..1; 0..1; 0..1] or read from file
  if(!params->read_grid_from_file)     
    CURVILINEAR_GRID<T>::Init_Local_Grid_Node_Positions_With_Uniform_Cube();
  else
    mpi_driver->Read_Global_Array_From_Disk(params->grid_filename, *grid);
  //CURVILINEAR_GRID<T>::Init_Local_Grid_Node_Positions_From_File();

  //mpi_driver->Write_Global_Array_To_Disk("node_position", *grid);

  //position moving grid at t = 0
  Set_Initial_Moving_Grid_Node_Positions();  

  //local grid: [xmin..xmax; ymin..ymax; zmin..zmax] + streching, if any
  CURVILINEAR_GRID<T>::Custom_Adjust_Local_Grid_Node_Positions();
  CURVILINEAR_GRID<T>::Finish_Initialization_Based_On_Node_Positions();

  Jacobian_diff = new ARRAY_3D<T>(*grid); 
  U_grid_xi = new ARRAY_3D<T>(*grid);
  U_grid_et = new ARRAY_3D<T>(*grid);
  U_grid_zt = new ARRAY_3D<T>(*grid);
  //grid copy from previous timestep
  prev_grid = new ARRAY_3D<VECTOR_3D<T> >(*grid);
  prev_XI = new ARRAY_3D<VECTOR_3D<T> >(*grid);
  prev_ET = new ARRAY_3D<VECTOR_3D<T> >(*grid);
  prev_ZT = new ARRAY_3D<VECTOR_3D<T> >(*grid);
  prev_U_grid_xi = new ARRAY_3D<T>(*grid);
  prev_U_grid_et = new ARRAY_3D<T>(*grid);
  prev_U_grid_zt = new ARRAY_3D<T>(*grid);

  //recording volume change with time
  volume_evolution = new T[parameters->max_timestep];
  vol_array_size = 0;
} 
//*****************************************************************************
// Destructor
//*****************************************************************************
template<class T>
CURVILINEAR_MOVING_GRID<T>::~CURVILINEAR_MOVING_GRID()
{
  delete Jacobian_diff; delete U_grid_xi; delete U_grid_et; delete U_grid_zt;
  delete prev_U_grid_xi; delete prev_U_grid_et; delete prev_U_grid_zt; 
  delete prev_grid; delete prev_XI; delete prev_ET; delete prev_ZT;
  delete[] volume_evolution;
}
//*****************************************************************************
// Calculate grid node velocities based on new node locations
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Update_Node_Velocities()
{
  // calculate new node velocities
  T coef = (T).25 / parameters->delta_time, 
    one_third = (T)1/(T)3, two_third = (T)2/(T)3;
  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo()-1; i++)
    for(int j=J_Min_With_Halo(); j<=J_Max_With_Halo()-1; j++)
      for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo()-1; k++){
	VECTOR_3D<T> dXdXI = (*grid)(i  ,j,k) - (*prev_grid)(i  ,j,k) 
	                   + (*grid)(i+1,j,k) - (*prev_grid)(i+1,j,k),
	             dXdET = (*grid)(i,j  ,k) - (*prev_grid)(i,j  ,k) 
	                   + (*grid)(i,j+1,k) - (*prev_grid)(i,j+1,k),
	             dXdZT = (*grid)(i,j,k  ) - (*prev_grid)(i,j,k  ) 
	                   + (*grid)(i,j,k+1) - (*prev_grid)(i,j,k+1),
	             XI((*XI_x)(i,j,k), (*XI_y)(i,j,k), (*XI_z)(i,j,k)), 
	             ET((*ET_x)(i,j,k), (*ET_y)(i,j,k), (*ET_z)(i,j,k)), 
	             ZT((*ZT_x)(i,j,k), (*ZT_y)(i,j,k), (*ZT_z)(i,j,k)),
	             XI_half_t(XI+(*prev_XI)(i,j,k)), 
                     ET_half_t(ET+(*prev_ET)(i,j,k)),
	             ZT_half_t(ZT+(*prev_ZT)(i,j,k));

	(*U_grid_xi)(i,j,k)= coef * VECTOR_3D<T>::Dot_Product(XI_half_t, dXdXI);
	(*U_grid_et)(i,j,k)= coef * VECTOR_3D<T>::Dot_Product(ET_half_t, dXdET);
	(*U_grid_zt)(i,j,k)= 0.;//coef * VECTOR_3D<T>::Dot_Product(ZT_half_t, dXdZT);
       	// use Adams-Bashforth for all timesteps after 1st
	if(parameters->time_step > 1){
	  (*U_grid_xi)(i,j,k) = two_third *(*U_grid_xi)(i,j,k)
                              + one_third *(*prev_U_grid_xi)(i,j,k);
	  (*U_grid_et)(i,j,k) = two_third *(*U_grid_et)(i,j,k)
                              + one_third *(*prev_U_grid_et)(i,j,k);
	  (*U_grid_zt)(i,j,k) = two_third *(*U_grid_zt)(i,j,k)
                              + one_third *(*prev_U_grid_zt)(i,j,k);
	}//else 
	  //(*U_grid_xi)(i,j,k)=(*U_grid_et)(i,j,k)=(*U_grid_zt)(i,j,k) = (T)0;	
  }

  /* alternative velocity calculation  
  ARRAY_3D<T> v_grid_old(*U_grid_et);
  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo(); i++)
    for(int j=J_Min_With_Halo()+1; j<=J_Max_With_Halo(); j++)
      for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo(); k++){
        (*U_grid_xi)(i,J_Min_With_Halo(),k) = 
        (*U_grid_zt)(i,J_Min_With_Halo(),k) = 
	(*U_grid_xi)(i,j,k) = (*U_grid_zt)(i,j,k) = (T)0;
	if(parameters->time_step > 1)
	  (*U_grid_et)(i,j,k) = (*U_grid_et)(i,j-1,k) 
                              + 2./3.*(*Jacobian_diff)(i,j,k)
	                      + 1./3.*(v_grid_old(i,j,k)-v_grid_old(i,j-1,k));
        else
   	  (*U_grid_et)(i,j,k) = (*U_grid_et)(i,j-1,k)+(*Jacobian_diff)(i,j,k);
      }
  */
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*U_grid_xi);
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*U_grid_et);
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*U_grid_zt);
}
//*****************************************************************************
// Calculate difference in cell volumes based on grid velocity fluxes
// Jacobian_diff: J^n+1 - J^n
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Update_Jacobian_Difference()
{
  for(int i=I_Min_With_Halo()+1; i<=I_Max_With_Halo(); i++)
    for(int j=J_Min_With_Halo()+1; j<=J_Max_With_Halo(); j++)
      for(int k=K_Min_With_Halo()+1; k<=K_Max_With_Halo(); k++){
	if(parameters->time_step > 1)
	  (*Jacobian_diff)(i,j,k) = //changed sign to -
	    - ((T)1.5*( (*U_grid_xi)(i-1,j,k) - (*U_grid_xi)(i,j,k) 
                      + (*U_grid_et)(i,j-1,k) - (*U_grid_et)(i,j,k)
                      + (*U_grid_zt)(i,j,k-1) - (*U_grid_zt)(i,j,k))
	      -(T)0.5*( (*prev_U_grid_xi)(i-1,j,k) - (*prev_U_grid_xi)(i,j,k) 
		      + (*prev_U_grid_et)(i,j-1,k) - (*prev_U_grid_et)(i,j,k)
		      + (*prev_U_grid_zt)(i,j,k-1) - (*prev_U_grid_zt)(i,j,k)));
        else
	  (*Jacobian_diff)(i,j,k) = -((*U_grid_xi)(i-1,j,k) -(*U_grid_xi)(i,j,k)
                                 + (*U_grid_et)(i,j-1,k) - (*U_grid_et)(i,j,k)
	                         + (*U_grid_zt)(i,j,k-1) - (*U_grid_zt)(i,j,k));
        /*added update of Jacobian*/
	(*inverse_Jacobian)(i,j,k) +=
                               parameters->delta_time * (*Jacobian_diff)(i,j,k);
      }
  //mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*Jacobian_diff);
}
//*****************************************************************************
// Calculate V of domain top & bottom boundaries based on new depth location
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Update_Bed_And_Lid_Velocities()
{
  assert(parameters->bed_velocity); assert(parameters->lid_velocity); 

  T velocity_sum = (T)0, 
    hz = CURVILINEAR_GRID<T>::Halo_Size(),
    i_size_w_h = CURVILINEAR_GRID<T>::I_Size()+2*hz, 
    k_size_w_h = CURVILINEAR_GRID<T>::K_Size()+2*hz;

  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo(); i++)
    for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo(); k++){
            (*parameters->bed_velocity)(i,k) = VECTOR_3D<T>((T)0,
                                      (  (*grid)(i,J_Min_With_Halo(),k).y 
                                       - (*prev_grid)(i,J_Min_With_Halo(),k).y )
                                      / parameters->delta_time,
//      - ((*parameters->depth)(i,k) - old_depth(i,k)) / parameters->delta_time,
						      (T)0);
      velocity_sum += (*parameters->bed_velocity)(i,k).y;
  }
  velocity_sum /= i_size_w_h*k_size_w_h; // calculate average bed velocity

  parameters->lid_velocity->Set_All_Elements_To(
  				       VECTOR_3D<T>((T)0, velocity_sum, (T)0));
}  
//*****************************************************************************
// Set node position based on the oscillating grid example (ex#1 from YJ paper)
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Set_Node_Position_Example_1(T current_time)
{ 
  //local domain min boundary
  int x_local_min = mpi_driver->my_coords_in_grid[0] * i_size,
      y_local_min = mpi_driver->my_coords_in_grid[1] * j_size, 
      z_local_min = mpi_driver->my_coords_in_grid[2] * k_size;

  T period = (T)400,
    t_over_T = current_time / (period * parameters->delta_time);

  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo(); i++)
    for(int j=J_Min_With_Halo(); j<=J_Max_With_Halo(); j++)
      for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo(); k++){
	T X_unif   = (T(x_local_min+i) - .5) / (T)parameters->num_total_nodes_x,
	  Y_unif   = (T(y_local_min+j) - .5) / (T)parameters->num_total_nodes_y,
	  Z_unif   = (T(z_local_min+k) - .5) / (T)parameters->num_total_nodes_z,
	  A_x = (T)1 + sin((T)2*parameters->pi*(Z_unif+t_over_T)),
	  A_y = (T)1 + sin((T)2*parameters->pi*(Z_unif+t_over_T)),

	  X_new =  x_length *(exp(A_x*X_unif) - (T)1) / (exp(A_x) - (T)1),
	  Y_new =  y_length * ((T)1 - (exp(A_y*Y_unif)-(T)1) / (exp(A_y)-(T)1)),
	  A_z = (T)1 + cos((T)2*parameters->pi*(X_new/x_length + t_over_T))
	             + cos((T)2*parameters->pi*(Y_new/y_length + t_over_T)),
	  Z_new = -z_length * ((T)1 - (exp(A_z*Z_unif)-(T)1) / (exp(A_z)-(T)1));
	(*grid)(i,j,k) = VECTOR_3D<T>(X_new, Y_new, Z_new);
  }
}
//*****************************************************************************
// Set grid node positions at t=0
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Set_Initial_Moving_Grid_Node_Positions()
{
  //Set_Node_Position_Example_1((T)0); //Example#1
}
//*****************************************************************************
// Adjust node positions to keep cells, being pulled apart, intact
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Fix_New_Grid()
{
  ARRAY_3D<T> vertical_grad(*grid);

  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo(); i++) 
    for(int j=J_Min_With_Halo(); j<=J_Max_With_Halo(); j++)
      for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo(); k++)
	vertical_grad(i,j,k) = (*grid)(i,j,k).y - (*prev_grid)(i,j,k).y;
  
  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo()-1; i++)
    for(int j=J_Min_With_Halo(); j<=J_Max_With_Halo()-1; j++)
      for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo()-1; k++)
	if((vertical_grad(i,j  ,k) < (T)0 || vertical_grad(i+1,j  ,k)< (T)0) && 
           (vertical_grad(i,j+1,k) > (T)0 || vertical_grad(i+1,j+1,k)> (T)0) ){
	  (*grid)(i  ,j  ,k) = (*prev_grid)(i  ,j  ,k);
	  (*grid)(i+1,j  ,k) = (*prev_grid)(i+1,j  ,k);
	  (*grid)(i  ,j+1,k) = (*prev_grid)(i  ,j+1,k);
	  (*grid)(i+1,j+1,k) = (*prev_grid)(i+1,j+1,k);
	}
}
//*****************************************************************************
// Record total domain volume evolution with time
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Write_Volume_Evolution_To_Disk()
{
  if(!mpi_driver->my_rank)
    mpi_driver->Write_Local_Array_To_Disk("volume_evolution", volume_evolution, 
					 vol_array_size, parameters->time_step);
}
//*****************************************************************************
// Save current grid before modifying it
//*****************************************************************************
template<class T>
void CURVILINEAR_MOVING_GRID<T>::Save_Grid_From_Previous_Timestep()
{
  // save domain volume
  volume_evolution[vol_array_size++] = 
                                 CURVILINEAR_GRID<T>::Calculate_Domain_Volume();
  // save current node positions, fluid and grid velocity fluxes, jacobians
  for(int i=I_Min_With_Halo(); i<=I_Max_With_Halo(); i++)
    for(int j=J_Min_With_Halo(); j<=J_Max_With_Halo(); j++)
      for(int k=K_Min_With_Halo(); k<=K_Max_With_Halo(); k++){
	(*prev_grid)(i,j,k) = (*grid)(i,j,k);
	(*prev_XI)(i,j,k) = 
                   VECTOR_3D<T>((*XI_x)(i,j,k), (*XI_y)(i,j,k), (*XI_z)(i,j,k));
	(*prev_ET)(i,j,k) = 
                   VECTOR_3D<T>((*ET_x)(i,j,k), (*ET_y)(i,j,k), (*ET_z)(i,j,k));
	(*prev_ZT)(i,j,k) = 
                   VECTOR_3D<T>((*ZT_x)(i,j,k), (*ZT_y)(i,j,k), (*ZT_z)(i,j,k));
	(*prev_U_grid_xi)(i,j,k) = (*U_grid_xi)(i,j,k);
	(*prev_U_grid_et)(i,j,k) = (*U_grid_et)(i,j,k);
	(*prev_U_grid_zt)(i,j,k) = (*U_grid_zt)(i,j,k);
      }
}
//*****************************************************************************
#endif
