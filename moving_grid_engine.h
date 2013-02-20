//*****************************************************************************
// File:  moving_grid_engine.h
// Class: MOVING_GRID_ENGINE<T>
// Note:  Instrumentation of moving grid logic
//*****************************************************************************
#ifndef __MOVING_GRID_ENGINE__
#define __MOVING_GRID_ENGINE__

#include "mpi_driver.h"
#include "interpolant.h"

template<class T=double>
class MOVING_GRID_ENGINE
{
 public:
   MOVING_GRID_ENGINE(PARAMETERS<T> *p, MPI_DRIVER<T> *mpi, CONVECTION<T> *c,
    CURVILINEAR_MOVING_GRID<T> *mg, ARRAY_3D<T> *r, ARRAY_3D<VECTOR_3D<T> > *u);
  ~MOVING_GRID_ENGINE();

  void Move_Grid_Nodes();

 private:
  void Update_Node_Velocities(ARRAY_3D<VECTOR_3D<T> >& old_grid, 
              ARRAY_3D<VECTOR_3D<T> >& old_XI, ARRAY_3D<VECTOR_3D<T> >& old_ET,
                                               ARRAY_3D<VECTOR_3D<T> >& old_ZT);
  void Update_Jacobian_Difference(ARRAY_3D<VECTOR_3D<T> >& old_grid, 
                        ARRAY_3D<T>& old_U_grid_xi, ARRAY_3D<T>& old_U_grid_et, 
				                    ARRAY_3D<T>& old_U_grid_zt);
  void Update_Bed_And_Lid_Velocities(ARRAY_2D<T>& old_depth);

  void Set_Node_Position_Example_1(T current_time);
  void Move_Nodes_Vertically_Based_On_Fluid_Velocity();

  void Evolve_Nodes_With_MMPDE();
  void Update_Monitor_Function();
  void Smooth_Monitor_Function_In_Space();
  void Smooth_Monitor_Function_In_Time();

  PARAMETERS<T> *parameters;
  MPI_DRIVER<T> *mpi_driver;
  CONVECTION<T> *convection;
  CURVILINEAR_MOVING_GRID<T> *moving_grid;
  ARRAY_3D<T> *rho, *M, *prev_M; // monitor function
  ARRAY_3D<VECTOR_3D<T> > *U_fluid, *prev_U_fluid;
  int imin_w_h, imax_w_h, jmin_w_h, jmax_w_h, kmin_w_h, kmax_w_h, halo;
  int number_of_smoothing_iters;
};
//*****************************************************************************
// Constructor
//*****************************************************************************
template<class T>
MOVING_GRID_ENGINE<T>::MOVING_GRID_ENGINE(PARAMETERS<T> *p, MPI_DRIVER<T> *mpi,
                             CONVECTION<T> *c, CURVILINEAR_MOVING_GRID<T> *mg,
                             ARRAY_3D<T> *r, ARRAY_3D<VECTOR_3D<T> > *u) 
 : parameters(p), mpi_driver(mpi), convection(c), 
   moving_grid(mg), rho(r), U_fluid(u)
{
  M = new ARRAY_3D<T>(*rho);  prev_M = NULL;
  number_of_smoothing_iters = 100;//(512le)//10;(256x64-le)//5;(128x32-le) //10;(sloshing)
  //copy of fluid velocity from previous time_step (for AB2) 
  prev_U_fluid = new ARRAY_3D<VECTOR_3D<T> >(*U_fluid);

  imin_w_h = moving_grid->I_Min_With_Halo(); 
  imax_w_h = moving_grid->I_Max_With_Halo();
  jmin_w_h = moving_grid->J_Min_With_Halo(); 
  jmax_w_h = moving_grid->J_Max_With_Halo();
  kmin_w_h = moving_grid->K_Min_With_Halo(); 
  kmax_w_h = moving_grid->K_Max_With_Halo(); halo = moving_grid->Halo_Size();
}
//*****************************************************************************
// Destructor
//*****************************************************************************
template<class T>
MOVING_GRID_ENGINE<T>::~MOVING_GRID_ENGINE() 
{
  delete M; if(prev_M) delete prev_M; delete prev_U_fluid;
}
//*****************************************************************************
// Moves grid nodes and updates auxiliary quantities
//*****************************************************************************
template<class T>
void MOVING_GRID_ENGINE<T>::Move_Grid_Nodes()
{
  moving_grid->Save_Grid_From_Previous_Timestep();
  //move grid nodes of new_grid based on some scheme
  //Move_Nodes_Based_On_Fluid_Velocity();
  //if(parameters->time_step % 150 == 0) // || parameters->time_step==1) 
  if(parameters->time_step>1)
  //    Evolve_Nodes_With_MMPDE();
  //else
    Move_Nodes_Vertically_Based_On_Fluid_Velocity();
  //mpi_driver->Write_Global_Array_To_Disk("node_position", *moving_grid->grid);
  moving_grid->Calculate_Metrics();
  moving_grid->Create_Update_Subgrids(false); // just update (don't re-create)

  moving_grid->Update_Node_Velocities();
  //*moving_grid->U_grid_et = *convection->Get_U_et();

  moving_grid->Update_Jacobian_Difference();
}
//*****************************************************************************
template<class T>
void MOVING_GRID_ENGINE<T>::Evolve_Nodes_With_MMPDE()
{
  Update_Monitor_Function();
  Smooth_Monitor_Function_In_Space();
  //Smooth_Monitor_Function_In_Time();

  T tolerance, 
    tolerance_threshold =.05*parameters->x_length/parameters->num_total_nodes_x;
  int iter_num = 0, gs_iter_num = 0, max_gs_iters = 3, max_iters = 1;
  int k = kmin_w_h;
  ARRAY_3D<VECTOR_3D<T> > saved_grid(*moving_grid->grid), 
                          grid(*moving_grid->grid),
                          prev_grid(*moving_grid->grid);
  //INTERPOLANT<T> interpolant(mpi_driver, *moving_grid->grid, *M);

  do{
    //interpolate M from saved_grid -> grid for iter_num>0
    if(iter_num){
      //interpolant.Interpolate_On_New_Grid(grid, *M);
      saved_grid = grid;
      gs_iter_num = 0;
    }
    do{ //GS iteration loop   
      prev_grid = grid;   
      T jmin_boundary = jmin_w_h+1, jmax_boundary = jmax_w_h-1,
	imin_boundary = imin_w_h+1, imax_boundary = imax_w_h-1;
      if(mpi_driver->suth_proc == MPI_PROC_NULL) jmin_boundary=jmin_w_h+3;
      if(mpi_driver->nrth_proc == MPI_PROC_NULL) jmax_boundary=jmax_w_h-3;
      if(mpi_driver->west_proc == MPI_PROC_NULL) imin_boundary=imin_w_h+3;//le
      if(mpi_driver->east_proc == MPI_PROC_NULL) imax_boundary=imax_w_h-3;//le
	
      // move internal nodes
      //      for(int i=imin_w_h+1; i<=imax_w_h-1; i++) //sloshing
      for(int i=imin_boundary; i<=imax_boundary; i++)
	for(int j=jmin_boundary; j<=jmax_boundary; j++)
	//for(int k=kmin_w_h+1; k<=kmax_w_h-1; k++)
	  {
	    T alpha_ip_j = (T).5*((*M)(i+1,j,k) + (*M)(i,j,k)),
	      alpha_im_j = (T).5*((*M)(i-1,j,k) + (*M)(i,j,k)),
	      beta_i_jp  = (T).5*((*M)(i,j+1,k) + (*M)(i,j,k)),
	      beta_i_jm  = (T).5*((*M)(i,j-1,k) + (*M)(i,j,k)),
	      denom = alpha_ip_j + alpha_im_j + beta_i_jp + beta_i_jm;
	    VECTOR_3D<T> rhs = alpha_ip_j * prev_grid(i+1,j,k) + 
	                       beta_i_jp  * prev_grid(i,j+1,k);
	    //sloshing: changed to .y only (!)
	    grid(i,j,k).x = ( rhs.x + alpha_im_j * grid(i-1,j,k).x
	    	                    + beta_i_jm  * grid(i,j-1,k).x ) / denom;
	    grid(i,j,k).y = ( rhs.y + alpha_im_j * grid(i-1,j,k).y
	 	                    + beta_i_jm  * grid(i,j-1,k).y ) / denom;
	}  
      // move boundary nodes
      // along South and North boundaries
      // changed to .x only (!)
      // changed jmin_w_h -> jmin_boundary-1, jmax_w_h -> jmax_boundary+1
      for(int i=imin_boundary; i<=imax_boundary; i++){
	grid(i,jmin_boundary-1,k).x += 
             grid(i,jmin_boundary,k).x - prev_grid(i,jmin_boundary,k).x; //suth
	grid(i,jmax_boundary+1,k).x +=
             grid(i,jmax_boundary,k).x - prev_grid(i,jmax_boundary,k).x; //nrth
	if(mpi_driver->suth_proc == MPI_PROC_NULL) grid(i,jmin_w_h,k).x = 
                                grid(i,jmin_w_h+1,k).x = grid(i,jmin_w_h+2,k).x;
	if(mpi_driver->nrth_proc == MPI_PROC_NULL) grid(i,jmax_w_h,k).x = 
			       grid(i,jmax_w_h-1,k).x =  grid(i,jmax_w_h-2,k).x;
      }
      
      // along West and East boundaries
      // changed to .y only (!)
      // changed imin_w_h -> imin_boundary-1, imax_w_h -> imax_boundary+1
      for(int j=jmin_boundary; j<=jmax_boundary; j++){
	grid(imin_boundary-1,j,k).y += 
	  grid(imin_boundary,j,k).y - prev_grid(imin_boundary,j,k).y; //west
	grid(imax_boundary+1,j,k).y += 
	  grid(imax_boundary,j,k).y - prev_grid(imax_boundary,j,k).y; //east
	if(mpi_driver->west_proc == MPI_PROC_NULL) grid(imin_w_h,j,k).y = 
		 	     grid(imin_w_h+1,j,k).y = grid(imin_w_h+2,j,k).y;
	if(mpi_driver->east_proc == MPI_PROC_NULL) grid(imax_w_h,j,k).y = 
			     grid(imax_w_h-1,j,k).y = grid(imax_w_h-2,j,k).y;

      }
      
      // back, front ...
 
      // calculate tolerance
      tolerance = (T)0;
      for(int i=imin_w_h; i<=imax_w_h; i++)
	for(int j=jmin_w_h; j<=jmax_w_h; j++)
	  //for(int k=kmin; k<=kmax; k++)
	  tolerance=fmax(tolerance, (grid(i,j,k)-prev_grid(i,j,k)).Magnitude());
      gs_iter_num++;
    }while(tolerance > tolerance_threshold && gs_iter_num < max_gs_iters);
    // calculate tolerance after all GS iterations
    tolerance = (T)0;
    for(int i=imin_w_h; i<=imax_w_h; i++)
      for(int j=jmin_w_h; j<=jmax_w_h; j++)
	//for(int k=kmin; k<=kmax; k++)
	tolerance =fmax(tolerance, (grid(i,j,k)-saved_grid(i,j,k)).Magnitude());
    iter_num++;
  }while(tolerance > tolerance_threshold && iter_num < max_iters);
  // do not change anything in case iteration count is == 1
  mpi_driver->Replace_With_Max_Value_Among_All_Procs(gs_iter_num),
  mpi_driver->Replace_With_Max_Value_Among_All_Procs(iter_num);
  if(!mpi_driver->my_rank) cout<<"GS iterations "<<gs_iter_num<<", "
                               <<"Outer iterations "<<iter_num<<endl;
  if(gs_iter_num>1 || iter_num>1)
    // replicate to all z-slices (do not change z coordinate)
    for(int i=imin_w_h; i<=imax_w_h; i++)
      for(int j=jmin_w_h; j<=jmax_w_h; j++)
	for(int kc=kmin_w_h; kc<=kmax_w_h; kc++){
	  (*moving_grid)(i,j,kc).x = grid(i,j,k).x;
	  (*moving_grid)(i,j,kc).y = grid(i,j,k).y;
	}  
  // update halo regions
  mpi_driver->Exchange_Ghost_Values_For_Vector_Field(*moving_grid->grid);
}
//*****************************************************************************
template<class T>
void MOVING_GRID_ENGINE<T>::Update_Monitor_Function()
{
  T alpha = 1000;//(512)//.0000005;(512le)//000001;(256le)//.0000025;(128le)//.0001;(slosh)
  for(int i=imin_w_h; i<=imax_w_h; i++)
    for(int j=jmin_w_h; j<=jmax_w_h; j++)
      for(int k=kmin_w_h; k<=kmax_w_h; k++){
	int i_minus = (i != imin_w_h ? i-1 : i), 
	    i_plus  = (i != imax_w_h ? i+1 : i),
            j_minus = (j != jmin_w_h ? j-1 : j), 
	    j_plus  = (j != jmax_w_h ? j+1 : j),
            k_minus = (k != kmin_w_h ? k-1 : k), 
	    k_plus  = (k != kmax_w_h ? k+1 : k);
		
	T dCdX = ((*rho)(i_plus,j,k) - (*rho)(i_minus,j,k)) / 
	         ((*moving_grid)(i_plus,j,k).x - (*moving_grid)(i_minus,j,k).x),
	  dCdY = ((*rho)(i,j_plus,k) - (*rho)(i,j_minus,k)) / 
	         ((*moving_grid)(i,j_plus,k).y - (*moving_grid)(i,j_minus,k).y),
	  dCdZ = 0.;
	   // ((*rho)(i,j,k_plus) - (*rho)(i,j,k_minus)) / 
	   // ((*moving_grid)(i,j,k_plus).z - (*moving_grid)(i,j,k_minus).z);
	   /*
	T dCdXI = (*rho)(i_plus,j,k) - (*rho)(i_minus,j,k),
	  dCdET = (*rho)(i,j_plus,k) - (*rho)(i,j_minus,k),
	  dCdX = dCdXI * (*moving_grid->XI_x)(i,j,k)
	       + dCdET * (*moving_grid->ET_x)(i,j,k),
	  dCdY = dCdXI * (*moving_grid->XI_y)(i,j,k)
	       + dCdET * (*moving_grid->ET_y)(i,j,k),
	  dCdZ = (T)0;
	   */
        (*M)(i,j,k) = sqrt( (T)1 + alpha*(dCdX*dCdX + dCdY*dCdY + dCdZ*dCdZ) );
      }
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*M);
}
//*****************************************************************************
template<class T>
void MOVING_GRID_ENGINE<T>::Smooth_Monitor_Function_In_Space()
{
  /*T cw = 1./8., w1 = 1./16., w2 = 1./32., w3 = 1./64.;
  for(int iters=0; iters<number_of_smoothing_iters; iters++)
    for(int i=imin_w_h+1; i<=imax_w_h-1; i++)
      for(int j=jmin_w_h+1; j<=jmax_w_h-1; j++)
	for(int k=kmin_w_h+1; k<=kmax_w_h-1; k++){
	  VECTOR_3D<T> li = M(i-1,j,k), ri = M(i+1,j,k), 
	               lj = M(i,j-1,k), rj = M(i,j+1,k),
	               lk = M(i,j,k-1), rk = M(i,j,k+1),
	               lilj = M(i-1,j-1,k), lirj = M(i-1,j+1,k), 
                       rilj = M(i+1,j-1,k), rirj = M(i+1,j+1,k),
	               lilk = M(i-1,j,k-1), lirk = M(i-1,j,k+1), 
                       rilk = M(i+1,j,k-1), rirk = M(i+1,j,k+1),
	               lklj = M(i,j-1,k-1), lkrj = M(i,j+1,k-1), 
                       rklj = M(i,j-1,k+1), rkrj = M(i,j+1,k+1),
	               liljlk = M(i-1,j-1,k-1), lirjlk = M(i-1,j+1,k-1),
	               riljlk = M(i+1,j-1,k-1), rirjlk = M(i+1,j+1,k-1),
	               liljrk = M(i-1,j-1,k+1), lirjrk = M(i-1,j+1,k+1),
   	               riljrk = M(i+1,j-1,k+1), rirjrk = M(i+1,j+1,k+1);

	  (*M)(i,j,k) = cw*(*M)(i,j,k) + w1*(li + ri + lj + rj + lk + rk) +
	    w2*(lilj+lirj+rilj+rirj+ lilk+lirk+rilk+rirk+ lklj+lkrj+rklj+rkrj) +
	    w3*(liljlk+lirjlk+riljlk+rirjlk + liljrk+lirjrk+riljrk+rirjrk);
   */
  for(int iters=0; iters<number_of_smoothing_iters; iters++)
    for(int i=imin_w_h; i<=imax_w_h; i++)
      for(int j=jmin_w_h; j<=jmax_w_h; j++)
	for(int k=kmin_w_h; k<=kmax_w_h; k++){
	  int num_terms = 0;
	  T M_terms = (T)0; 
	  if(i>imin_w_h){
	    num_terms++;
	    M_terms += (*M)(i-1,j,k);
	  }
	  if(i<imax_w_h){
	    num_terms++;
	    M_terms += (*M)(i+1,j,k);
	  }
	  if(j>jmin_w_h){
	    num_terms++;
	    M_terms += (*M)(i,j-1,k);
	  }
	  if(j<jmax_w_h){
	    num_terms++;
	    M_terms += (*M)(i,j+1,k);
	  }
	  T inv_num_terms = 1./(T)num_terms; //.25 if 4 terms
	  (*M)(i,j,k) += inv_num_terms*M_terms;
	  (*M)(i,j,k) /= 2.;
	}
  // fill-in ghost cells based on updated boundary adjacent nodes
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*M);
}
//*****************************************************************************
template<class T>
void MOVING_GRID_ENGINE<T>::Smooth_Monitor_Function_In_Time()
{
  if(prev_M){
    *M = (T).8 * (*M) + (T).2 * (*prev_M);
    *prev_M = *M;
  } else prev_M = new ARRAY_3D<T>(*M);
}
//*****************************************************************************
// Move grid nodes based on AB2 scheme using fluid velocity
//*****************************************************************************
template<class T>
void MOVING_GRID_ENGINE<T>::Move_Nodes_Vertically_Based_On_Fluid_Velocity()
{
  int i_min = moving_grid->I_Min(), i_max = moving_grid->I_Max(), 
      j_min = moving_grid->J_Min(), j_max = moving_grid->J_Max(), 
      k_min = moving_grid->K_Min(), k_max = moving_grid->K_Max();
  if(mpi_driver->west_proc == MPI_PROC_NULL) i_min -= halo;
  if(mpi_driver->east_proc == MPI_PROC_NULL) i_max += halo;
  if(mpi_driver->suth_proc == MPI_PROC_NULL) j_min += 1; // don't move bottom
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) j_max -= 1; // don't move top
  if(mpi_driver->back_proc == MPI_PROC_NULL) k_min -= halo;
  if(mpi_driver->frnt_proc == MPI_PROC_NULL) k_max += halo;

  // if bed/lid are moving, move bottom/top ghost nodes
  if(parameters->bed_velocity && mpi_driver->suth_proc == MPI_PROC_NULL) 
    j_min = moving_grid->J_Min_With_Halo();
  if(parameters->lid_velocity && mpi_driver->nrth_proc == MPI_PROC_NULL) 
    j_max = moving_grid->J_Max_With_Halo();

  ARRAY_3D<T> grid_velocity(*U_fluid,2);
  //smooth velocity in y
  T weight = 1; 
  int num_iters = 5;
  for(int it=1; it<=num_iters; it++)
  for(int i=moving_grid->I_Min_With_Halo(); i<=moving_grid->I_Max_With_Halo(); i++)
    for(int j=moving_grid->J_Min(); j<=moving_grid->J_Max(); j++)
      for(int k=moving_grid->K_Min_With_Halo(); k<=moving_grid->K_Max_With_Halo(); k++){
	int terms = 4;
	/*if(j>moving_grid->J_Min()){
	  terms++;
	  grid_velocity(i,j,k) += weight * grid_velocity(i,j-3,k);
	}
	if(j<moving_grid->J_Max()){
	  terms++;
	  grid_velocity(i,j,k) += weight * grid_velocity(i,j+3,k);
	  }*/
	grid_velocity(i,j,k) += weight * (
 	                       grid_velocity(i,j-2,k) + grid_velocity(i,j-1,k) 
 			     + grid_velocity(i,j+1,k) + grid_velocity(i,j+2,k));

        grid_velocity(i,j,k) /= 1. + (T)terms * weight;
      }
  //smooth velocity in x
  for(int it=1; it<=num_iters; it++)
  for(int i=moving_grid->I_Min(); i<=moving_grid->I_Max(); i++)
    for(int j=moving_grid->J_Min_With_Halo(); j<=moving_grid->J_Max_With_Halo(); j++)
      for(int k=moving_grid->K_Min_With_Halo(); k<=moving_grid->K_Max_With_Halo(); k++){
	int terms = 4;
	/*
	if(i>moving_grid->I_Min()){
	  terms++;
	  grid_velocity(i,j,k) += weight * grid_velocity(i-3,j,k);
	}
	if(i<moving_grid->I_Max()){
	  terms++;
	  grid_velocity(i,j,k) += weight * grid_velocity(i+3,j,k);
	}
	*/
	grid_velocity(i,j,k) += weight * (
	                       grid_velocity(i-2,j,k) + grid_velocity(i-1,j,k) +
                             + grid_velocity(i+1,j,k) + grid_velocity(i+2,j,k));
        grid_velocity(i,j,k) /= 1. + (T)terms * weight;
      }
  mpi_driver->Write_Global_Array_To_Disk("smooth_fluid_v", grid_velocity,
					  parameters->time_step);
  /*
  for(int i=moving_grid->I_Min(); i<=moving_grid->I_Max(); i++)
    for(int j=moving_grid->J_Min(); j<=moving_grid->J_Max(); j++)
      for(int k=moving_grid->K_Min(); k<=moving_grid->K_Max(); k++){
	(*prev_U_fluid)(i,j,k).y += (*prev_U_fluid)(i-1,j-1,k).y 
                                  + (*prev_U_fluid)(i+1,j-1,k).y 
                                  + (*prev_U_fluid)(i-1,j+1,k).y 
                                  + (*prev_U_fluid)(i+1,j+1,k).y;
	(*prev_U_fluid)(i,j,k).y /= 5.;
      }
 */

  // move internal nodes
  for(int i=i_min; i<=i_max; i++)
    for(int j=j_min; j<=j_max; j++)
      for(int k=k_min; k<=k_max; k++)
	if(parameters->time_step>1){
	 	  (*moving_grid)(i,j,k).y += parameters->delta_time*
                                  ( (T)1.5*grid_velocity(i,j,k)
				    //( (T)1.5*(*U_fluid)(i,j,k).y 
				    - (T).5*(*prev_U_fluid)(i,j,k).y );
	  /*moving_grid)(i,j,k) += parameters->delta_time*
                                  ( (T)1.5*(*U_fluid)(i,j,k) 
                                  - (T).5*(*prev_U_fluid)(i,j,k) );*/
          (*prev_U_fluid)(i,j,k).y = grid_velocity(i,j,k);
	}else
  	  (*moving_grid)(i,j,k).y += parameters->delta_time*(*U_fluid)(i,j,k).y;
  /* move boundaries: West - East
  i_min = moving_grid->I_Min(); i_max = moving_grid->I_Max();
  j_min = moving_grid->J_Min(); j_max = moving_grid->J_Max();
  k_min = moving_grid->K_Min(); k_max = moving_grid->K_Max();
  for(int j=j_min-2; j<=j_max+2; j++)
    for(int k=k_min-2; k<=k_max+2; k++)
      for(int i=1; i<=2; i++){
	(*moving_grid)(i_min-i,j,k).y += 
	   (*moving_grid)(i_min,j,k).y - (*moving_grid->prev_grid)(i_min,j,k).y;
	(*moving_grid)(i_max+i,j,k).y += 
	   (*moving_grid)(i_max,j,k).y - (*moving_grid->prev_grid)(i_max,j,k).y;
      }
  */
  // fill-in ghost cells based on updated boundary adjacent nodes
  mpi_driver->Exchange_Ghost_Values_For_Vector_Field(*moving_grid->grid);
 
  // save current fluid velocity
  //*prev_U_fluid = *U_fluid;
}
//*****************************************************************************
#endif
