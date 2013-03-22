// convection.h
// Performs momentum and scalar convection

#ifndef __CONVECTION__
#define __CONVECTION__

#include "curvilinear_grid.h"
#include "universal_limiter.h"

template<class T=double>
class CONVECTION
{
 public:
  CONVECTION(PARAMETERS<T> *p, MPI_DRIVER<T> *md, CURVILINEAR_GRID<T> *g, 
            ARRAY_3D<T> *rhoin, ARRAY_3D<VECTOR_3D<T> > *u,
            ARRAY_3D<T> *uxi, ARRAY_3D<T> *uet, ARRAY_3D<T> *uzt) 
 : grid(g), parameters(p), mpi_driver(md), 
   U(u), U_xi(uxi), U_et(uet), U_zt(uzt), rho(rhoin), 
   lower_boundary(md->local_grid_lower_bound), 
   upper_boundary(md->local_grid_upper_bound)
  {  
    if(parameters->moving_grid){
      U_grid_xi = ((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_xi;
      U_grid_et = ((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_et;
      U_grid_zt = ((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_zt;
    }else 
      U_grid_xi = U_grid_et = U_grid_et = NULL;  

    switch(parameters->universal_limiter){
      case UPWIND: limiter = new UPWIND_LIMITER<T>();break;
      case MUSCL: limiter = new MUSCL_LIMITER<T>();break;
      default: limiter = NULL;break;
    }
  }
  ~CONVECTION(){delete limiter;}

  void Add_Quick_Scheme_Convection_Term(ARRAY_3D<VECTOR_3D<T> >& expression);
  void Add_Scalar_Convection_Term(ARRAY_3D<T>& expression);
 
  void Quick_Velocity_Flux_Update(ARRAY_3D<VECTOR_3D<T> >& u);
  void Central_Velocity_Flux_Update(ARRAY_3D<VECTOR_3D<T> >& u);
  void Sharp(ARRAY_3D<T>& scalar, ARRAY_3D<VECTOR_3D<T> >& scalar_flux,
                      ARRAY_3D<T>& U_xi, ARRAY_3D<T>& U_et, ARRAY_3D<T>& U_zt);
  void TVD(ARRAY_3D<T>& scalar, ARRAY_3D<VECTOR_3D<T> >& scalar_flux,
                      ARRAY_3D<T>& U_xi, ARRAY_3D<T>& U_et, ARRAY_3D<T>& U_zt);
 
  void Add_Moving_Grid_Convection_Term(ARRAY_3D<VECTOR_3D<T> >& expression);
  void Add_Moving_Grid_Convection_Term(ARRAY_3D<T>& expression);

  ARRAY_3D<T>* Get_U_xi() { return U_xi;}
  ARRAY_3D<T>* Get_U_et() { return U_et;}
  ARRAY_3D<T>* Get_U_zt() { return U_zt;}

  CURVILINEAR_GRID<T> *grid;
 private:
  PARAMETERS<T> *parameters;

  UNIVERSAL_LIMITER<T> *limiter;
  MPI_DRIVER<T> *mpi_driver;
  ARRAY_3D<VECTOR_3D<T> > *U;
  ARRAY_3D<T> *rho;
  ARRAY_3D<T> *U_xi, *U_et, *U_zt,                // fluxes on faces  
              *U_grid_xi, *U_grid_et, *U_grid_zt; 
  int *lower_boundary, *upper_boundary;           // set per proc
};
//*****************************************************************************
// Adds a convective term based on the QUICK scheme to 'expression'
//*****************************************************************************
template<class T>
void CONVECTION<T>::Add_Quick_Scheme_Convection_Term(
				            ARRAY_3D<VECTOR_3D<T> >& expression)
{
  //assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(*grid->grid,expression));
  T f_neg_east, f_pos_east, f_neg_west, f_pos_west, 
    f_neg_nrth, f_pos_nrth, f_neg_suth, f_pos_suth,
    f_neg_frnt, f_pos_frnt, f_neg_back, f_pos_back, rms,
    a_east_u, a_west_u, a_nrth_u, a_suth_u, a_frnt_u, a_back_u, a_cntr_u;

  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
	f_neg_east = (T).5* ((*U_xi)(i,  j,  k) - fabs((*U_xi)(i,  j,  k  )) );
        f_pos_east = (T).5* ((*U_xi)(i,  j,  k) + fabs((*U_xi)(i,  j,  k  )) );
	f_neg_west = (T).5* ((*U_xi)(i-1,j,  k) - fabs((*U_xi)(i-1,j,  k  )) );
	f_pos_west = (T).5* ((*U_xi)(i-1,j,  k) + fabs((*U_xi)(i-1,j,  k  )) );
	f_neg_nrth = (T).5* ((*U_et)(i,  j,  k) - fabs((*U_et)(i,  j,  k  )) );
	f_pos_nrth = (T).5* ((*U_et)(i,  j,  k) + fabs((*U_et)(i,  j,  k  )) );
	f_neg_suth = (T).5* ((*U_et)(i,  j-1,k) - fabs((*U_et)(i,  j-1,k  )) );
	f_pos_suth = (T).5* ((*U_et)(i,  j-1,k) + fabs((*U_et)(i,  j-1,k  )) );
	f_neg_frnt = (T).5* ((*U_zt)(i,  j,  k) - fabs((*U_zt)(i,  j,  k  )) );
	f_pos_frnt = (T).5* ((*U_zt)(i,  j,  k) + fabs((*U_zt)(i,  j,  k  )) );
	f_neg_back = (T).5* ((*U_zt)(i,  j,k-1) - fabs((*U_zt)(i,  j,  k-1)) );
	f_pos_back = (T).5* ((*U_zt)(i,  j,k-1) + fabs((*U_zt)(i,  j,  k-1)) );

	rms = (*U_xi)(i,j,k)-(*U_xi)(i-1,j,k) + (*U_et)(i,j,k)-(*U_et)(i,j-1,k)
            + (*U_zt)(i,j,k)-(*U_zt)(i,j,k-1);

	a_east_u  = (T)-0.5 * (*U_xi)(i,j,k) 
	          + (T).125 * (f_pos_east - (T)2*f_neg_east - f_neg_west);
	a_west_u  = (T).5   * (*U_xi)(i-1,j,k) 
	          + (T).125 * (f_pos_east + (T)2*f_pos_west - f_neg_west);
	a_nrth_u  = (T)-0.5 * (*U_et)(i,j,k)
	          + (T).125 * (f_pos_nrth - (T)2*f_neg_nrth - f_neg_suth);
	a_suth_u  = (T).5   * (*U_et)(i,j-1,k)
	          + (T).125 * (f_pos_nrth + (T)2*f_pos_suth - f_neg_suth);
	a_frnt_u  = (T)-.5  * (*U_zt)(i,j,k)
	          + (T).125 * (f_pos_frnt - (T)2*f_neg_frnt - f_neg_back);
	a_back_u  = (T).5   * (*U_zt)(i,j,k-1)
	          + (T).125 * (f_pos_frnt + (T)2*f_pos_back - f_neg_back );
	a_cntr_u  = 
            a_east_u + a_west_u + a_nrth_u + a_suth_u + a_frnt_u + a_back_u
	  + (T).125 * ( f_neg_east + f_neg_nrth + f_neg_frnt 
                      - f_pos_west - f_pos_suth - f_pos_back  ) + rms;
	// add convective term
	expression(i,j,k) +=a_east_u * (*U)(i+1,j,k) + a_west_u * (*U)(i-1,j,k) 
	                  + a_nrth_u * (*U)(i,j+1,k) + a_suth_u * (*U)(i,j-1,k)
 	                  + a_frnt_u * (*U)(i,j,k+1) + a_back_u * (*U)(i,j,k-1)
	                  - a_cntr_u * (*U)(i,j,k)
       + (T).125 * (f_neg_east * (*U)(i+2,j,k) - f_pos_west  * (*U)(i-2,j,k)
        	  + f_neg_nrth * (*U)(i,j+2,k) - f_pos_suth  * (*U)(i,j-2,k)
		  + f_neg_frnt * (*U)(i,j,k+2) - f_pos_back  * (*U)(i,j,k-2));
      }
} 
//*****************************************************************************
// Adds a convective term based on the SHARP or TVD scheme to 'expression'
//*****************************************************************************
template<class T>
void CONVECTION<T>::Add_Scalar_Convection_Term(ARRAY_3D<T>& expression)
{
  //assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(*grid->grid, expression));

  ARRAY_3D<VECTOR_3D<T> > rho_flux(*rho);
  if(parameters->universal_limiter == SHARP) 
    Sharp(*rho, rho_flux, *U_xi, *U_et, *U_zt);
  else
    TVD(*rho, rho_flux, *U_xi, *U_et, *U_zt);
  //mpi_driver->Write_Global_Array_To_Disk("rho_flux", rho_flux, 
  //					 parameters->time_step);
  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++)
	expression(i,j,k) += 
	   (*rho)(i,j,k) * ( (*U_xi)(i, j, k) - (*U_xi)(i-1,j  ,k  )
			     + (*U_et)(i, j, k) - (*U_et)(i  ,j-1,k  )	
		    	     + (*U_zt)(i, j, k) - (*U_zt)(i  ,j  ,k-1) )

	                   - ( (*U_xi)(i  , j  ,k  )*rho_flux(i  ,j  ,k  ).x 
                             - (*U_xi)(i-1, j  ,k  )*rho_flux(i-1,j  ,k  ).x
	                     + (*U_et)(i  , j  ,k  )*rho_flux(i  ,j  ,k  ).y 
	                     - (*U_et)(i  , j-1,k  )*rho_flux(i  ,j-1,k  ).y 
	                     + (*U_zt)(i  , j  ,k  )*rho_flux(i  ,j  ,k  ).z 
	                     - (*U_zt)(i  , j  ,k-1)*rho_flux(i  ,j  ,k-1).z );
}
//*****************************************************************************
// QUICK interpolation: updates values of velocity fluxes U_xi, U_et, U_zt
//*****************************************************************************
template<class T>
void CONVECTION<T>::Quick_Velocity_Flux_Update(ARRAY_3D<VECTOR_3D<T> >& u)
{  
  assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(u, *grid->grid));
  // Calculating U_xi
  for(int i = lower_boundary[0]; i <= upper_boundary[0]; i++) 
    //for(int i = U_xi->I_Min()-1; i <= U_xi->I_Max(); i++)
    for(int j = U_xi->J_Min(); j <= U_xi->J_Max(); j++)
      for(int k = U_xi->K_Min(); k <= U_xi->K_Max(); k++)
        if((*U_xi)(i,j,k) >= (T)0)
          (*U_xi)(i,j,k) = (*grid->XI_x)(i,j,k) *
            ((T).75*u(i,j,k).x + (T).375*u(i+1,j,k).x - (T).125*u(i-1,j,k).x)
            + (*grid->XI_y)(i,j,k) *
            ((T).75*u(i,j,k).y + (T).375*u(i+1,j,k).y - (T).125*u(i-1,j,k).y)
            + (*grid->XI_z)(i,j,k) *
            ((T).75*u(i,j,k).z + (T).375*u(i+1,j,k).z - (T).125*u(i-1,j,k).z);
        else
          (*U_xi)(i,j,k) = (*grid->XI_x)(i,j,k) *
            ((T).75*u(i+1,j,k).x + (T).375*u(i,j,k).x - (T).125*u(i+2,j,k).x)
            + (*grid->XI_y)(i,j,k) *
            ((T).75*u(i+1,j,k).y + (T).375*u(i,j,k).y - (T).125*u(i+2,j,k).y)
            + (*grid->XI_z)(i,j,k) *
            ((T).75*u(i+1,j,k).z + (T).375*u(i,j,k).z - (T).125*u(i+2,j,k).z);

  // Calculating U_et
  for(int i = U_et->I_Min(); i <= U_et->I_Max(); i++) 
    for(int j = lower_boundary[1]; j <= upper_boundary[1]; j++)
      //for(int j = U_et->J_Min()-1; j <= U_et->J_Max(); j++)
      for(int k = U_et->K_Min(); k <= U_et->K_Max(); k++)
        if((*U_et)(i,j,k) >= (T)0)
          (*U_et)(i,j,k) = (*grid->ET_x)(i,j,k) *
            ((T).75*u(i,j,k).x + (T).375*u(i,j+1,k).x - (T).125*u(i,j-1,k).x)
            + (*grid->ET_y)(i,j,k) *
            ((T).75*u(i,j,k).y + (T).375*u(i,j+1,k).y - (T).125*u(i,j-1,k).y)
            + (*grid->ET_z)(i,j,k)*
            ((T).75*u(i,j,k).z + (T).375*u(i,j+1,k).z - (T).125*u(i,j-1,k).z);
        else
          (*U_et)(i,j,k) = (*grid->ET_x)(i,j,k) *
            ((T).75*u(i,j+1,k).x + (T).375*u(i,j,k).x - (T).125*u(i,j+2,k).x)
            + (*grid->ET_y)(i,j,k) *
            ((T).75*u(i,j+1,k).y + (T).375*u(i,j,k).y - (T).125*u(i,j+2,k).y)
            + (*grid->ET_z)(i,j,k) *
            ((T).75*u(i,j+1,k).z + (T).375*u(i,j,k).z - (T).125*u(i,j+2,k).z);

  // Update flux on upper boundary based on lid_velocity
  if(parameters->lid_velocity)
    for(int i = U_et->I_Min(); i <= U_et->I_Max(); i++)
      for(int k = U_et->K_Min(); k <= U_et->K_Max(); k++)
        (*U_et)(i,upper_boundary[1]+1,k) = 
          (*grid->ET_y)(i,upper_boundary[1]+1,k)
          * (*parameters->lid_velocity)(i,k).y;

  // Update flux on bottom boundary based on bed_velocity
  if(parameters->bed_velocity)
    for(int i = U_et->I_Min(); i <= U_et->I_Max(); i++)
      for(int k = U_et->K_Min(); k <= U_et->K_Max(); k++)
        (*U_et)(i,lower_boundary[1]-1,k) = 
          (*grid->ET_y)(i,lower_boundary[1]-1,k)
          * (*parameters->bed_velocity)(i,k).y;

  // Update flux on west boundary based on west_velocity
  T Qp = 0., Q = 0., alpha; 

  if(parameters->west_velocity){
    if(mpi_driver->west_proc == MPI_PROC_NULL){
      for(int j = U_xi->J_Min(); j <= U_xi->J_Max(); j++)
        for(int k = U_xi->K_Min(); k <= U_xi->K_Max(); k++){
          (*U_xi)(lower_boundary[0]-1,j,k) = 
            (*grid->XI_x)(lower_boundary[0]-1,j,k)
            * (*parameters->west_velocity)(j,k).x;
          Qp += (*U_xi)(lower_boundary[0]-1,j,k);
        }
    }
  }

  //Adjust to conserve volume
  mpi_driver->Replace_With_Sum_On_All_Procs(Qp);

  if(parameters->west_velocity){
    if(mpi_driver->west_proc == MPI_PROC_NULL){
      alpha = -Qp/parameters->num_total_nodes_y/parameters->num_total_nodes_z;
      for(int j = U_xi->J_Min(); j <= U_xi->J_Max(); j++)
        for(int k = U_xi->K_Min(); k <= U_xi->K_Max(); k++){
          (*U_xi)(lower_boundary[0]-1,j,k) += alpha;
          Q += (*U_xi)(lower_boundary[0]-1,j,k);
        }
    }
  }

  // Calculating U_zt
  for(int i = U_zt->I_Min(); i <= U_zt->I_Max(); i++) 
    for(int j = U_zt->J_Min(); j <= U_zt->J_Max(); j++)
      //for(int k = U_zt->K_Min()-1; k <= U_zt->K_Max(); k++)
      for(int k = lower_boundary[2]; k <= upper_boundary[2]; k++)
        if((*U_zt)(i,j,k) >= (T)0)
          (*U_zt)(i,j,k) = (*grid->ZT_x)(i,j,k) *
            ((T).75*u(i,j,k).x + (T).375*u(i,j,k+1).x - (T).125*u(i,j,k-1).x)
            + (*grid->ZT_y)(i,j,k) *
            ((T).75*u(i,j,k).y + (T).375*u(i,j,k+1).y - (T).125*u(i,j,k-1).y)
            + (*grid->ZT_z)(i,j,k) *
            ((T).75*u(i,j,k).z + (T).375*u(i,j,k+1).z - (T).125*u(i,j,k-1).z);
        else
          (*U_zt)(i,j,k) = (*grid->ZT_x)(i,j,k) *
            ((T).75*u(i,j,k+1).x + (T).375*u(i,j,k).x - (T).125*u(i,j,k+2).x)
            + (*grid->ZT_y)(i,j,k) *
            ((T).75*u(i,j,k+1).y + (T).375*u(i,j,k).y - (T).125*u(i,j,k+2).y)
            + (*grid->ZT_z)(i,j,k) *
            ((T).75*u(i,j,k+1).z + (T).375*u(i,j,k).z - (T).125*u(i,j,k+2).z);
}
//*****************************************************************************
// Central Difference interpolation: updates velocity fluxes U_xi, U_et, U_zt
//*****************************************************************************
template<class T>
void CONVECTION<T>::Central_Velocity_Flux_Update(ARRAY_3D<VECTOR_3D<T> >& u)
{  
  assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(u, *grid->grid));
  // Calculating U_xi
  for(int i = U_xi->I_Min()-1; i <= U_xi->I_Max(); i++) 
    for(int j = U_xi->J_Min(); j <= U_xi->J_Max(); j++)
      for(int k = U_xi->K_Min(); k <= U_xi->K_Max(); k++)
	(*U_xi)(i,j,k) = .5*( (*grid->XI_x)(i,j,k)*(u(i+1,j,k).x+u(i,j,k).x)
	                    + (*grid->XI_y)(i,j,k)*(u(i+1,j,k).y+u(i,j,k).y)
			    + (*grid->XI_z)(i,j,k)*(u(i+1,j,k).z+u(i,j,k).z));
  // Calculating U_et
  for(int i = U_et->I_Min(); i <= U_et->I_Max(); i++) 
    for(int j = U_xi->J_Min()-1; j <= U_xi->J_Max(); j++)
      for(int k = U_et->K_Min(); k <= U_et->K_Max(); k++)
	(*U_et)(i,j,k) = .5*( (*grid->ET_x)(i,j,k)*(u(i,j+1,k).x+u(i,j,k).x)
	                    + (*grid->ET_y)(i,j,k)*(u(i,j+1,k).y+u(i,j,k).y)
			    + (*grid->ET_z)(i,j,k)*(u(i,j+1,k).z+u(i,j,k).z));

  // Update flux on upper boundary based on lid_velocity
  if(parameters->lid_velocity)
    for(int i = U_et->I_Min(); i <= U_et->I_Max(); i++)
      for(int k = U_et->K_Min(); k <= U_et->K_Max(); k++)
	(*U_et)(i,upper_boundary[1]+1,k) = 
	                               (*grid->ET_y)(i,upper_boundary[1]+1,k)
	                             * (*parameters->lid_velocity)(i,k).y;

  // Update flux on bottom boundary based on bed_velocity
  if(parameters->bed_velocity)
    for(int i = U_et->I_Min(); i <= U_et->I_Max(); i++)
      for(int k = U_et->K_Min(); k <= U_et->K_Max(); k++)
	(*U_et)(i,lower_boundary[1]-1,k) = 
	                               (*grid->ET_y)(i,lower_boundary[1]-1,k)
	                             * (*parameters->bed_velocity)(i,k).y;
  // Calculating U_zt
  for(int i = U_zt->I_Min(); i <= U_zt->I_Max(); i++) 
    for(int j = U_zt->J_Min(); j <= U_zt->J_Max(); j++)
      for(int k = U_xi->K_Min()-1; k <= U_xi->K_Max(); k++)
	(*U_zt)(i,j,k) = .5*( (*grid->ZT_x)(i,j,k)*(u(i,j,k+1).x+u(i,j,k).x)
	                    + (*grid->ZT_y)(i,j,k)*(u(i,j,k+1).y+u(i,j,k).y)
			    + (*grid->ZT_z)(i,j,k)*(u(i,j,k+1).z+u(i,j,k).z));
}
//*****************************************************************************
// SHARP advection scheme: 
// calculates scalar_fluxes based on scalar and face velocities (U_xi,U_et,U_zt)
//*****************************************************************************
template<class T>
void CONVECTION<T>::Sharp(ARRAY_3D<T>& scalar, 
			  ARRAY_3D<VECTOR_3D<T> >& scalar_flux,
                        ARRAY_3D<T>& U_xi, ARRAY_3D<T>& U_et, ARRAY_3D<T>& U_zt)
{  
  assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(*grid->grid, scalar));
  assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(scalar_flux, scalar));
  int x_min=scalar.I_Min(), x_max=scalar.I_Max(),
      y_min=scalar.J_Min(), y_max=scalar.J_Max(),
      z_min=scalar.K_Min(), z_max=scalar.K_Max();
  T tolerance = 1e-5;
  ARRAY_3D<VECTOR_3D<T> > cf(x_min,x_max,y_min,y_max,z_min,z_max,1);

  //calculate I-direction
  for(int i = x_min-1; i <= x_max; i++)
    for(int j = y_min-1; j <= y_max; j++)
      for(int k = z_min-1; k <= z_max; k++)
	if(U_xi(i,j,k) >= (T)0) {
	  T delf = scalar(i+1,j,k) - scalar(i-1,j,k);
	  if(abs(delf) < tolerance) cf(i,j,k).x = (T)0.125;
	  else{ 
	    T scalar_east = (scalar(i,j,k) - scalar(i-1,j,k)) / delf;
	    if(scalar_east <= (T)-1 || scalar_east >= (T)1.5) 
	      cf(i,j,k).x = (T)0.125;
	    else
	      if(scalar_east <= (T).25)
		if(scalar_east <= (T)0) 
		  cf(i,j,k).x = (T).5 + (T).375 * scalar_east;
		else cf(i,j,k).x = (T).5 - (T).625 * sqrt(scalar_east);
	      else
		if(scalar_east <= (T)1) 
		  cf(i,j,k).x = (T).25 * ((T)1 - scalar_east);
		else cf(i,j,k).x = (T)-.25 * ((T)1 - scalar_east);
	  }
	  scalar_flux(i,j,k).x = (T).5 * (scalar(i,j,k)+scalar(i+1,j,k))
             - cf(i,j,k).x * (scalar(i-1,j,k)-2.*scalar(i,j,k)+scalar(i+1,j,k));
	}else{
	  T delf = scalar(i,j,k) - scalar(i+2,j,k);
	  if( abs(delf) < tolerance) cf(i,j,k).x = (T)0.125;
	  else{ 
	    T scalar_east = (scalar(i+1,j,k) - scalar(i+2,j,k)) / delf;
	    if(scalar_east <= (T)-1 || scalar_east >= (T)1.5) 
	      cf(i,j,k).x = (T)0.125;
	    else
	      if(scalar_east <= (T).25)
		if(scalar_east <= (T)0) 
		  cf(i,j,k).x = (T).5 + (T).375 * scalar_east;
		else cf(i,j,k).x = (T).5 - (T).625 * sqrt(scalar_east);
	      else
		if(scalar_east <= (T)1) 
		  cf(i,j,k).x = (T).25 * ((T)1 - scalar_east);
		else cf(i,j,k).x = (T)-.25 * ((T)1 - scalar_east);
	  }
	  scalar_flux(i,j,k).x = (T).5 * (scalar(i,j,k)+scalar(i+1,j,k))
	      - cf(i,j,k).x * (scalar(i,j,k)-2*scalar(i+1,j,k)+scalar(i+2,j,k));
	}

  //calculate J-direction
  for(int i = x_min-1; i <= x_max; i++)
    for(int j = y_min-1; j <= y_max; j++)
      for(int k = z_min-1; k <= z_max; k++)
	if(U_et(i,j,k) >= (T)0) {
	  T delta = scalar(i,j+1,k) - scalar(i,j-1,k);
	  if(abs(delta) < tolerance) cf(i,j,k).y = (T)0.125;
	  else{ 
	    T scalar_east = (scalar(i,j,k) - scalar(i,j-1,k)) / delta;
	    if(scalar_east <= (T)-1 || scalar_east >= (T)1.5) 
	      cf(i,j,k).y = (T)0.125;
	    else
	      if(scalar_east <= (T).25)
		if(scalar_east <= (T)0) 
		  cf(i,j,k).y = (T).5 + (T).375 * scalar_east;
		else 
		  cf(i,j,k).y = (T).5 - (T).625 * sqrt(scalar_east);
	      else
		if(scalar_east <= (T)1) 
		  cf(i,j,k).y =  (T).25 * ((T)1 - scalar_east);
		else 
		  cf(i,j,k).y = (T)-.25 * ((T)1 - scalar_east);
	  }
	  scalar_flux(i,j,k).y = (T).5 * (scalar(i,j,k)+scalar(i,j+1,k)) 
	      - cf(i,j,k).y * (scalar(i,j-1,k)-2*scalar(i,j,k)+scalar(i,j+1,k));
	}else{
	  T delta = scalar(i,j,k) - scalar(i,j+2,k);
	  if(abs(delta) < tolerance) cf(i,j,k).y = (T)0.125;
	  else{ 
	    T scalar_east = (scalar(i,j+1,k) - scalar(i,j+2,k)) / delta;
	    if(scalar_east <= (T)-1 || scalar_east >= (T)1.5) 
	      cf(i,j,k).y = (T)0.125;
	    else
	      if(scalar_east <= (T).25)
		if(scalar_east <= (T)0) 
		  cf(i,j,k).y = (T).5 + (T).375 * scalar_east;
		else 
		  cf(i,j,k).y = (T).5 - (T).625 * sqrt(scalar_east);
	      else
		if(scalar_east <= (T)1) 
		  cf(i,j,k).y =  (T).25 * ((T)1 - scalar_east);
		else 
		  cf(i,j,k).y = (T)-.25 * ((T)1 - scalar_east);
	  }
	  scalar_flux(i,j,k).y = (T).5 * (scalar(i,j,k) + scalar(i,j+1,k))
	     - cf(i,j,k).y * (scalar(i,j,k)-2*scalar(i,j+1,k)+scalar(i,j+2,k));
	}

  //calculate K-direction
  for(int i = x_min-1; i <= x_max; i++)
    for(int j = y_min-1; j <= y_max; j++)
      for(int k = z_min-1; k <= z_max; k++)
	if(U_zt(i,j,k) >= (T)0) {
	  T delta = scalar(i,j,k+1) - scalar(i,j,k-1);
	  if(abs(delta) < tolerance) cf(i,j,k).z = (T)0.125;
	  else{ 
	    T scalar_east = (scalar(i,j,k) - scalar(i,j,k-1)) / delta;
	    if(scalar_east <= (T)-1 || scalar_east >= (T)1.5) 
	      cf(i,j,k).z = (T)0.125;
	    else
	      if(scalar_east <= (T).25)
		if(scalar_east <= (T)0) 
		  cf(i,j,k).z = (T).5 + (T).375 * scalar_east;
		else 
		  cf(i,j,k).z = (T).5 - (T).625 * sqrt(scalar_east);
	      else
		if(scalar_east <= (T)1) 
		  cf(i,j,k).z =  (T).25 * ((T)1 - scalar_east);
		else 
		  cf(i,j,k).z = (T)-.25 * ((T)1 - scalar_east);
	  }
	  scalar_flux(i,j,k).z = (T).5 * (scalar(i,j,k)+scalar(i,j,k+1))
	      - cf(i,j,k).z * (scalar(i,j,k-1)-2*scalar(i,j,k)+scalar(i,j,k+1));
	}else{
	  T delta = scalar(i,j,k) - scalar(i,j,k+2);
	  if(abs(delta) < tolerance) cf(i,j,k).z = (T)0.125;
	  else{ 
	    T scalar_east = (scalar(i,j,k+1) - scalar(i,j,k+2)) / delta;
	    if(scalar_east <= (T)-1 || scalar_east >= (T)1.5) 
	      cf(i,j,k).z = (T)0.125;
	    else
	      if(scalar_east <= (T).25)
		if(scalar_east <= (T)0) 
		  cf(i,j,k).z = (T).5 + (T).375 * scalar_east;
		else 
		  cf(i,j,k).z = (T).5 - (T).625 * sqrt(scalar_east);
	      else
		if(scalar_east <= (T)1) 
		  cf(i,j,k).z =  (T).25 * ((T)1 - scalar_east);
		else 
		  cf(i,j,k).z = (T)-.25 * ((T)1 - scalar_east);
	  }
	  scalar_flux(i,j,k).z = (T).5 * (scalar(i,j,k)+scalar(i,j,k+1)) 
	     - cf(i,j,k).z * (scalar(i,j,k)-2*scalar(i,j,k+1)+scalar(i,j,k+2));
	}	
}
//*****************************************************************************
// TVD advection scheme: 
// calculates scalar_fluxes based on scalar and face velocities (U_xi,U_et,U_zt)
// uses the limiter base class that is polymorphic based on the input parameter
//*****************************************************************************
template<class T>
void CONVECTION<T>::TVD(ARRAY_3D<T>& scalar, 
			  ARRAY_3D<VECTOR_3D<T> >& scalar_flux,
                        ARRAY_3D<T>& U_xi, ARRAY_3D<T>& U_et, ARRAY_3D<T>& U_zt)
{  
  assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(*grid->grid, scalar));
  assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(scalar_flux, scalar));
  int x_min=scalar.I_Min(), x_max=scalar.I_Max(),
      y_min=scalar.J_Min(), y_max=scalar.J_Max(),
      z_min=scalar.K_Min(), z_max=scalar.K_Max(), halo = scalar.Halo_Size();

  //calculate I-direction
  for(int i = x_min-1; i <= x_max; i++)
    for(int j = y_min-1; j <= y_max; j++)
      for(int k = z_min-1; k <= z_max; k++)
	if(U_xi(i,j,k) >= (T)0) {
	  T r_plus = (scalar(i  ,j,k) - scalar(i-1,j,k)) /
	             (scalar(i+1,j,k) - scalar(i  ,j,k)),	   
	    U_face = U_xi(i,j,k);//.5 * ((*U)(i,j,k)+(*U)(i+1,j,k));	  
	  scalar_flux(i,j,k).x = scalar(i,j,k) + (T).5*(*limiter)(r_plus)*
                                     (1-U_face)*(scalar(i+1,j,k)-scalar(i,j,k));
	}else{
	  T r_minus = (scalar(i+2,j,k) - scalar(i+1,j,k)) /
	              (scalar(i+1,j,k) - scalar(i  ,j,k)),
	    U_face = U_xi(i,j,k);//.5 * ((*U)(i,j,k)+(*U)(i+1,j,k));	  
	  scalar_flux(i,j,k).x = scalar(i+1,j,k) - (T).5*(*limiter)(r_minus)*
                                     (1+U_face)*(scalar(i+1,j,k)-scalar(i,j,k));
	}

  //calculate J-direction
  for(int i = x_min-1; i <= x_max; i++)
    for(int j = y_min-1; j <= y_max; j++)
      for(int k = z_min-1; k <= z_max; k++)
	if(U_et(i,j,k) >= (T)0) {
	  T r_plus = (scalar(i,j  ,k) - scalar(i,j-1,k)) /
	             (scalar(i,j+1,k) - scalar(i,j  ,k)),
	    U_face = U_et(i,j,k);//.5 * ((*U)(i,j,k)+(*U)(i,j+1,k));
	  scalar_flux(i,j,k).y = scalar(i,j,k) + (T).5*(*limiter)(r_plus)*
                                     (1-U_face)*(scalar(i,j+1,k)-scalar(i,j,k));
	}else{
	  T r_minus = (scalar(i,j+2,k) - scalar(i,j+1,k)) /
	              (scalar(i,j+1,k) - scalar(i,j  ,k)),	   
	    U_face = U_et(i,j,k);//.5 * ((*U)(i,j,k)+(*U)(i,j+1,k));	  
	  scalar_flux(i,j,k).y = scalar(i,j+1,k) - (T).5*(*limiter)(r_minus)*
                                     (1+U_face)*(scalar(i,j+1,k)-scalar(i,j,k));
	}

  //calculate K-direction
  for(int i = x_min-1; i <= x_max; i++)
    for(int j = y_min-1; j <= y_max; j++)
      for(int k = z_min-1; k <= z_max; k++)
	if(U_zt(i,j,k) >= (T)0) {
	  T r_plus = (scalar(i,j,k  ) - scalar(i,j,k-1)) /
	             (scalar(i,j,k+1) - scalar(i,j,k  )),
	    U_face = U_zt(i,j,k);//.5 * ((*U)(i,j,k)+(*U)(i,j,k+1));	  
	  scalar_flux(i,j,k).z = scalar(i,j,k) + (T).5*(*limiter)(r_plus)*
                                     (1-U_face)*(scalar(i,j,k+1)-scalar(i,j,k));
	}else{
	  T r_minus = (scalar(i,j,k+2) - scalar(i,j,k+1)) /
	              (scalar(i,j,k+1) - scalar(i,j,k  )),
	    U_face = U_zt(i,j,k);//.5 * ((*U)(i,j,k)+(*U)(i,j,k+1));	  
	  scalar_flux(i,j,k).z = scalar(i,j,k+1) - (T).5*(*limiter)(r_minus)*
                                     (1+U_face)*(scalar(i,j,k+1)-scalar(i,j,k));
	}	
}
//*****************************************************************************
// Adds a convective term for moving grid to VECTOR_3D 'expression'
// Added to the RHS: +div(u_i*U_g)
//*****************************************************************************
template<class T>
void CONVECTION<T>::Add_Moving_Grid_Convection_Term(
					    ARRAY_3D<VECTOR_3D<T> >& expression)
{
  //assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(expression,*grid->grid));
  assert(parameters->moving_grid);

  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++)
	expression(i,j,k) += 
	    (*U_grid_xi)(i  ,j  ,k  ) * .5 * ((*U)(i,j,k)+(*U)(i+1,j,k))
	  - (*U_grid_xi)(i-1,j  ,k  ) * .5 * ((*U)(i,j,k)+(*U)(i-1,j,k))
          + (*U_grid_et)(i  ,j  ,k  ) * .5 * ((*U)(i,j,k)+(*U)(i,j+1,k))
	  - (*U_grid_et)(i  ,j-1,k  ) * .5 * ((*U)(i,j,k)+(*U)(i,j-1,k))
          + (*U_grid_zt)(i  ,j  ,k  ) * .5 * ((*U)(i,j,k)+(*U)(i,j,k+1)) 
          - (*U_grid_zt)(i  ,j  ,k-1) * .5 * ((*U)(i,j,k)+(*U)(i,j,k-1)); 
}
//*****************************************************************************
// Adds(positive on RHS) a convective term for moving grid to scalar'expression'
// On RHS this term is: +div(rho*U_g)
//*****************************************************************************
template<class T>
void CONVECTION<T>::Add_Moving_Grid_Convection_Term(ARRAY_3D<T>& expression)
{
  //assert(ARRAY_3D<VECTOR_3D<T> >::Equal_Dimensions(*grid->grid, expression));
  assert(parameters->moving_grid);

  ARRAY_3D<VECTOR_3D<T> > rho_flux(*rho); 
  ARRAY_3D<T> nUg_xi((T)-1*(*U_grid_xi)), 
              nUg_et((T)-1*(*U_grid_et)),
              nUg_zt((T)-1*(*U_grid_zt));

   if(parameters->universal_limiter == SHARP) 
     Sharp(*rho, rho_flux, *U_grid_xi, *U_grid_et, *U_grid_zt);
     //Sharp(*rho, rho_flux, nUg_xi, nUg_et, nUg_zt);
   else
     TVD(*rho, rho_flux, *U_grid_xi, *U_grid_et, *U_grid_zt);
     //TVD(*rho, rho_flux, nUg_xi, nUg_et, nUg_zt);

  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++)
	expression(i,j,k) += 
 	  	   (*U_grid_xi)(i  , j  ,k  )*rho_flux(i  ,j  ,k  ).x 
                 - (*U_grid_xi)(i-1, j  ,k  )*rho_flux(i-1,j  ,k  ).x
	         + (*U_grid_et)(i  , j  ,k  )*rho_flux(i  ,j  ,k  ).y 
                 - (*U_grid_et)(i  , j-1,k  )*rho_flux(i  ,j-1,k  ).y 
                 + (*U_grid_zt)(i  , j  ,k  )*rho_flux(i  ,j  ,k  ).z 
	         - (*U_grid_zt)(i  , j  ,k-1)*rho_flux(i  ,j  ,k-1).z;
 /*
         (*U_grid_xi)(i  ,j  ,k  ) * .5 * ((*rho)(i,j,k)+(*rho)(i+1,j,k))
       - (*U_grid_xi)(i-1,j  ,k  ) * .5 * ((*rho)(i,j,k)+(*rho)(i-1,j,k))
       + (*U_grid_et)(i  ,j  ,k  ) * .5 * ((*rho)(i,j,k)+(*rho)(i,j+1,k))
       - (*U_grid_et)(i  ,j-1,k  ) * .5 * ((*rho)(i,j,k)+(*rho)(i,j-1,k))
       + (*U_grid_zt)(i  ,j  ,k  ) * .5 * ((*rho)(i,j,k)+(*rho)(i,j,k+1))
       - (*U_grid_zt)(i  ,j  ,k-1) * .5 * ((*rho)(i,j,k)+(*rho)(i,j,k-1));
 */  
}
//*****************************************************************************
#endif
