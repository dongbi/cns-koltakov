//*****************************************************************************
// Pressure solver class
//*****************************************************************************
#ifndef __PRESSURE__
#define __PRESSURE__

#include "curvilinear_grid.h"
#include "tridiagonal_solver.h"
#include "metric_quantities.h"
#include <cmath>

template<class T=double>
class PRESSURE
{
  public:
    PRESSURE(PARAMETERS<T> *par, MPI_DRIVER<T> *md, CURVILINEAR_GRID<T> *g, 
        ARRAY_3D<T> *p, ARRAY_3D<VECTOR_3D<T> > *vU, 
        ARRAY_3D<T> *uxi, ARRAY_3D<T> *uet, ARRAY_3D<T> *uzt);
    ~PRESSURE()
    {delete RHS; delete Residual;
      if(parameters->mg_sub_levels) {delete P_sub;delete RHS_sub;delete Res_sub;}
    }
    void Solve();

  private:
    void Reset_All_Arrays();
    void Condense_Array(ARRAY_3D<T>& P, ARRAY_3D<T>& P_condensed);
    void Extend_Array(ARRAY_3D<T>& P_condensed, ARRAY_3D<T>& P);
    void Compute_Residual(T& residual_l2, T& rhs_l2, T& residual_min,
        T& residual_max, ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Residual,
        METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new);
    void Smooth_Pressure(int level, T& residual_l2, T& rhs_l2, T& residual_min, 
        T& residual_max, ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Residual);
    void Smooth_Pressure_In_X(
        ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Res,
        METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new);  
    void Smooth_Pressure_In_Y(
        ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Res,
        METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new);  
    void Smooth_Pressure_In_Z(
        ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Res,
        METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new);
    void Fill_In_Residual_On_Six_Boundary_Faces(ARRAY_3D<T>& Res, ARRAY_3D<T>& P, 
        ARRAY_3D<T>& RHS, METRIC_QUANTITIES<T>& mq);  
    void Fill_In_Twelve_Boundary_Edges(ARRAY_3D<T>& Res);
    void Fill_In_Eight_Boundary_Corners(ARRAY_3D<T>& P);
    void Initialize_Metric_Quantities(int level, 
        METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new);
    CURVILINEAR_GRID<T> *grid;
    PARAMETERS<T> *parameters;
    MPI_DRIVER<T> *mpi_driver;

    ARRAY_3D<T> *P, *RHS, *Residual;
    ARRAY_3D<VECTOR_3D<T> > *u;
    ARRAY_3D<T> *U_xi, *U_et, *U_zt; 
    ARRAY_1D<ARRAY_3D<T>* > *P_sub, *RHS_sub, *Res_sub;
};
//*****************************************************************************
// Ctor
//*****************************************************************************
  template<class T>
PRESSURE<T>::PRESSURE(PARAMETERS<T> *par, MPI_DRIVER<T> *md, 
    CURVILINEAR_GRID<T> *g, ARRAY_3D<T> *p, ARRAY_3D<VECTOR_3D<T> > *vU, 
    ARRAY_3D<T> *uxi, ARRAY_3D<T> *uet, ARRAY_3D<T> *uzt)
: parameters(par), mpi_driver(md), grid(g), 
  P(p), u(vU), U_xi(uxi), U_et(uet), U_zt(uzt) 
{
  RHS = new ARRAY_3D<T>(*P,false);
  Residual = new ARRAY_3D<T>(*P,false);

  // create subgrid arrays in case there are any multigrid sublevels
  if(parameters->mg_sub_levels) {
    P_sub = new ARRAY_1D<ARRAY_3D<T>* >(parameters->mg_sub_levels);
    RHS_sub = new ARRAY_1D<ARRAY_3D<T>* >(parameters->mg_sub_levels);
    Res_sub = new ARRAY_1D<ARRAY_3D<T>* >(parameters->mg_sub_levels);
    for(int level = 1; level <=parameters->mg_sub_levels; level++) {      
      (*P_sub)(level)   = new ARRAY_3D<T>( (*grid->G11_sub)(level)->I_Min(),
          (*grid->G11_sub)(level)->I_Max(),
          (*grid->G11_sub)(level)->J_Min(), 
          (*grid->G11_sub)(level)->J_Max(), 
          (*grid->G11_sub)(level)->K_Min(), 
          (*grid->G11_sub)(level)->K_Max(), 
          P->Halo_Size() );
      (*RHS_sub)(level) = new ARRAY_3D<T>( *(*P_sub)(level), false);
      (*Res_sub)(level) = new ARRAY_3D<T>( *(*P_sub)(level), false);
    }
  }else
    P_sub = RHS_sub = Res_sub = NULL;
}  
//*****************************************************************************
// Helper function: sets all arrays to 0
//*****************************************************************************
template<class T> void PRESSURE<T>::Reset_All_Arrays()
{
  P->Set_All_Elements_To((T)0);
  RHS->Set_All_Elements_To((T)0);
  Residual->Set_All_Elements_To((T)0);
  for(int level = 1; level <=parameters->mg_sub_levels; level++) { 
    (*P_sub)(level)->Set_All_Elements_To((T)0);
    (*RHS_sub)(level)->Set_All_Elements_To((T)0);
    (*Res_sub)(level)->Set_All_Elements_To((T)0);
  }
}
//*****************************************************************************
// Main Poisson pressure solver function
//*****************************************************************************
template<class T> void PRESSURE<T>::Solve()
{
  T inv_dt = (T)1/parameters->delta_time,
    residual_l2 = (T)0, rhs_l2 = (T)0, residual_min = (T)0, residual_max = (T)0;

  Reset_All_Arrays();
  // Set RHS on 8 boundaries
  // RHS: west
  if(mpi_driver->west_proc == MPI_PROC_NULL)
    for(int j = RHS->J_Min(); j <= RHS->J_Max(); j++)
      for(int k = RHS->K_Min(); k <= RHS->K_Max(); k++) {
        int imin = RHS->I_Min();
        VECTOR_3D<T> v1 = VECTOR_3D<T>((*grid->XI_x)(imin-1,j,k),
            (*grid->XI_y)(imin-1,j,k),
            (*grid->XI_z)(imin-1,j,k)),
          v2 = (T)15*(*u)(imin,j,k)-(T)10*(*u)(imin+1,j,k)+(T)3*(*u)(imin+2,j,k);
        (*RHS)(imin-1,j,k) = (T).125*inv_dt * VECTOR_3D<T>::Dot_Product(v1,v2);
        /*
        //Account for inflow at western boundary
        if(parameters->west_velocity){
          if(mpi_driver->west_proc == MPI_PROC_NULL) {
            (*U_xi)(imin-1,j,k) = (*grid->XI_x)(imin-1,j,k)
              * (*parameters->west_velocity)(j,k).x;
          }
        }
        else {
          (*U_xi)(imin-1,j,k) = (T)0;
        }
        */
        if(!parameters->west_velocity) 
          (*U_xi)(imin-1,j,k) = (T)0;
      }
  // RHS: east
  if(mpi_driver->east_proc == MPI_PROC_NULL)
    for(int j = RHS->J_Min(); j <= RHS->J_Max(); j++)
      for(int k = RHS->K_Min(); k <= RHS->K_Max(); k++) {
        int imax = RHS->I_Max();
        VECTOR_3D<T> v1 = VECTOR_3D<T>( (*grid->XI_x)(imax,j,k),
            (*grid->XI_y)(imax,j,k),
            (*grid->XI_z)(imax,j,k) ),
          v2 = (T)15*(*u)(imax,j,k)-(T)10*(*u)(imax-1,j,k)+(T)3*(*u)(imax-2,j,k);
        (*RHS)(imax+1,j,k) = (T).125*inv_dt * VECTOR_3D<T>::Dot_Product(v1,v2);
        (*U_xi)(imax,j,k) = (T)0;
      }
  // RHS: front
  if(mpi_driver->frnt_proc == MPI_PROC_NULL)
    for(int i = RHS->I_Min(); i <= RHS->I_Max(); i++) 
      for(int k = RHS->K_Min(); k <= RHS->K_Max(); k++) {
        int jmin = RHS->J_Min();
        VECTOR_3D<T> v1 = VECTOR_3D<T>( (*grid->ET_x)(i,jmin-1,k),
            (*grid->ET_y)(i,jmin-1,k),
            (*grid->ET_z)(i,jmin-1,k) ),
          v2 = (T)15*(*u)(i,jmin,k)-(T)10*(*u)(i,jmin+1,k)+(T)3*(*u)(i,jmin+2,k);
        (*RHS)(i,jmin-1,k) = (T).125*inv_dt * VECTOR_3D<T>::Dot_Product(v1,v2);
        (*U_et)(i,jmin-1,k) = (T)0;
      }
  // RHS: back
  if(mpi_driver->back_proc == MPI_PROC_NULL)
    for(int i = RHS->I_Min(); i <= RHS->I_Max(); i++) 
      for(int k = RHS->K_Min(); k <= RHS->K_Max(); k++) {
        int jmax = RHS->J_Max();
        VECTOR_3D<T> v1 = VECTOR_3D<T>( (*grid->ET_x)(i,jmax,k),
            (*grid->ET_y)(i,jmax,k),
            (*grid->ET_z)(i,jmax,k) ),
          v2 = (T)15*(*u)(i,jmax,k)-(T)10*(*u)(i,jmax-1,k)+(T)3*(*u)(i,jmax-2,k);
        (*RHS)(i,jmax+1,k) = (T).125*inv_dt * VECTOR_3D<T>::Dot_Product(v1,v2);
        (*U_et)(i,jmax,k) = (T)0;
      }
  // RHS: south
  if(mpi_driver->suth_proc == MPI_PROC_NULL)
    for(int i = RHS->I_Min(); i <= RHS->I_Max(); i++) 
      for(int j = RHS->J_Min(); j <= RHS->J_Max(); j++) {
        int kmin = RHS->K_Min();
        VECTOR_3D<T> v1 = VECTOR_3D<T>( (*grid->ZT_x)(i,j,kmin-1),
            (*grid->ZT_y)(i,j,kmin-1),
            (*grid->ZT_z)(i,j,kmin-1) ),
          v2 = (T)15*(*u)(i,j,kmin)-(T)10*(*u)(i,j,kmin+1)+(T)3*(*u)(i,j,kmin+2);
        (*RHS)(i,j,kmin-1) = (T).125*inv_dt * VECTOR_3D<T>::Dot_Product(v1,v2);
        (*U_zt)(i,j,kmin-1) = (T)0;
      }
  // RHS: north
  if(mpi_driver->nrth_proc == MPI_PROC_NULL)
    for(int i = RHS->I_Min(); i <= RHS->I_Max(); i++) 
      for(int j = RHS->J_Min(); j <= RHS->J_Max(); j++) {
        int kmax = RHS->K_Max();
        VECTOR_3D<T> v1 = VECTOR_3D<T>( (*grid->ZT_x)(i,j,kmax),
            (*grid->ZT_y)(i,j,kmax),
            (*grid->ZT_z)(i,j,kmax) ),
          v2 = (T)15*(*u)(i,j,kmax)-(T)10*(*u)(i,j,kmax-1)+(T)3*(*u)(i,j,kmax-2);
        (*RHS)(i,j,kmax+1) = (T).125*inv_dt * VECTOR_3D<T>::Dot_Product(v1,v2);
        (*U_zt)(i,j,kmax) = (T)0;
      }
  // RHS: in the interior
  for(int i = RHS->I_Min(); i <= RHS->I_Max(); i++)
    for(int j = RHS->J_Min(); j <= RHS->J_Max(); j++)
      for(int k = RHS->K_Min(); k <= RHS->K_Max(); k++)
        (*RHS)(i,j,k) = inv_dt * ( (*U_xi)(i,j,k) - (*U_xi)(i-1,j,k)
            + (*U_et)(i,j,k) - (*U_et)(i,j-1,k)
            + (*U_zt)(i,j,k) - (*U_zt)(i,j,k-1) );
  // Main Iteration Loop  
  for(int iter = 1; iter <= parameters->max_mg_iters; iter++){

    // Coarsening (going down the V-cycle) if multigrid_levels!=0
    for(int level = 1; level <= parameters->mg_sub_levels; level++){
      if(level == 1){
        Smooth_Pressure(0, residual_l2,rhs_l2,residual_min,residual_max,*P,*RHS,
            *Residual); //initialize Residual 
        Condense_Array(*Residual, *(*RHS_sub)(1));
      } else
        Condense_Array(*(*Res_sub)(level-1), *(*RHS_sub)(level));
        (*P_sub)(level)->Set_All_Elements_To((T)0);
        Smooth_Pressure(level, residual_l2, rhs_l2, residual_min, residual_max, 
            *(*P_sub)(level), *(*RHS_sub)(level), *(*Res_sub)(level));
    }
    // Refinement (going up the V-cycle) if multigrid_levels!=0
    for(int level = parameters->mg_sub_levels; level >= 1; level--){
      if(level > 1){
        Extend_Array(*(*P_sub)(level), *(*P_sub)(level-1));
        Smooth_Pressure(level-1, residual_l2, rhs_l2, residual_min,residual_max,
            *(*P_sub)(level-1), *(*RHS_sub)(level-1), *(*Res_sub)(level-1));
      }
      else // last step (level==1)
        Extend_Array(*(*P_sub)(1), *P); //extend to finest level of P
    }
    // In case there are no MG sub-levels or last step of refinement loop 
    Smooth_Pressure(0, residual_l2, rhs_l2, residual_min, residual_max, 
        *P, *RHS, *Residual);
    // Check convergence
    T max_abs_error = fmax(fabs(residual_min), fabs(residual_max)),
      relative_resid = residual_l2 / (rhs_l2 < 1e-15 ? 1e-15 : rhs_l2);
    if(mpi_driver->my_rank==0  && iter%10==1 &&
        parameters->time_step % parameters->print_timestep_period==0) 
      cout<<"Iter = "<<iter
        <<", Rel_res = "<<relative_resid 
        <<", Rhs_l2 = "<<rhs_l2
        <<", Res_l2 = "<<residual_l2
        /*<<", Max_err = "<<max_abs_error*/<<endl;
    if( residual_l2    < parameters->mg_tol_absolute_resid &&
        max_abs_error  < parameters->mg_tol_error_resid    &&
        relative_resid < parameters->mg_tol_relative_resid ) {
      if(parameters->time_step % parameters->print_timestep_period==0 &&
          mpi_driver->my_rank==0) {
        cout<< "Total V-cycles of MG = " << iter << endl;
        cout<< "---" << endl;
      }
      break; //end iteration for-loop
    }
  }//for: iter
  //if(mpi_driver->my_rank==0) cout<<"Completed all "<<parameters->max_mg_iters
  //                               <<" V-cycles of MG" <<endl;
}
//*****************************************************************************
// Compresses P into P_condensed (half the size): fine to coarse
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Condense_Array(ARRAY_3D<T>& P, ARRAY_3D<T>& P_condensed)
{
  assert(P.Halo_Size()); assert(P_condensed.Halo_Size()); // halo > 0

  if(!parameters->two_d){
    // P_condensed: interior region
    for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
      for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic - 1, j = 2*jc - 1, k = 2*kc - 1; 
          P_condensed(ic,jc,kc) = P(i,j  ,k  ) + P(i+1,j  ,k  ) 
            + P(i,j+1,k  ) + P(i+1,j+1,k  )
            + P(i,j  ,k+1) + P(i+1,j  ,k+1)
            + P(i,j+1,k+1) + P(i+1,j+1,k+1);
        }
    if(mpi_driver->west_proc == MPI_PROC_NULL)
      for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = P.I_Min()-1, j = 2*jc - 1, k = 2*kc - 1; //i=0
          P_condensed(P_condensed.I_Min()-1, jc, kc) = P(i,j  ,k  ) 
            + P(i,j+1,k  )
            + P(i,j  ,k+1)
            + P(i,j+1,k+1);
        }	
    if(mpi_driver->east_proc == MPI_PROC_NULL)
      for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = P.I_Max()+1, j = 2*jc - 1, k = 2*kc - 1; //i=imax+1
          P_condensed(P_condensed.I_Max()+1, jc, kc) = P(i,j  ,k  ) 
            + P(i,j+1,k  )
            + P(i,j  ,k+1)
            + P(i,j+1,k+1);
        }	
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic - 1, j = P.J_Min()-1, k = 2*kc - 1; //j=0
          P_condensed(ic, P_condensed.J_Min()-1, kc) = P(i  ,j,k  ) 
            + P(i+1,j,k  )
            + P(i  ,j,k+1)
            + P(i+1,j,k+1);
        }	
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic - 1, j = P.J_Max()+1, k = 2*kc - 1; //j=jmax+1
          P_condensed(ic, P_condensed.J_Max()+1, kc) = P(i  ,j,k  ) 
            + P(i+1,j,k  )
            + P(i  ,j,k+1)
            + P(i+1,j,k+1);
        }
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++){
          int i = 2*ic - 1, j = 2*jc - 1, k = P.K_Min()-1; //k=0
          P_condensed(ic, jc, P_condensed.K_Min()-1) = P(i  ,j  ,k) 
            + P(i+1,j  ,k)
            + P(i  ,j+1,k)
            + P(i+1,j+1,k);
        }
    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++){
          int i = 2*ic - 1, j = 2*jc - 1, k = P.K_Max()+1; //k=kmax+1
          P_condensed(ic, jc, P_condensed.K_Max()+1) = P(i  ,j  ,k) 
            + P(i+1,j  ,k)
            + P(i  ,j+1,k)
            + P(i+1,j+1,k);
        } 
  } else {
    // P_condensed: interior region
    for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
      for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic - 1, j = jc, k = 2*kc - 1; 
          P_condensed(ic,jc,kc) = P(i,j  ,k  ) + P(i+1,j  ,k  ) 
            + P(i,j  ,k+1) + P(i+1,j  ,k+1);
        }
    if(mpi_driver->west_proc == MPI_PROC_NULL)
      for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = P.I_Min()-1, j = jc, k = 2*kc - 1; //i=0
          P_condensed(P_condensed.I_Min()-1, jc, kc) = P(i,j  ,k  ) 
            + P(i,j  ,k+1);
        }	
    if(mpi_driver->east_proc == MPI_PROC_NULL)
      for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = P.I_Max()+1, j = jc, k = 2*kc - 1; //i=imax+1
          P_condensed(P_condensed.I_Max()+1, jc, kc) = P(i,j  ,k  ) 
            + P(i,j  ,k+1);
        }	
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic - 1, j = P.J_Min()-1, k = 2*kc - 1; //j=0
          P_condensed(ic, P_condensed.J_Min()-1, kc) = P(i  ,j,k  ) 
            + P(i+1,j,k  )
            + P(i  ,j,k+1)
            + P(i+1,j,k+1);
        }	
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int kc = P_condensed.K_Min(); kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic - 1, j = P.J_Max()+1, k = 2*kc - 1; //j=jmax+1
          P_condensed(ic, P_condensed.J_Max()+1, kc) = P(i  ,j,k  ) 
            + P(i+1,j,k  )
            + P(i  ,j,k+1)
            + P(i+1,j,k+1);
        }
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++){
          int i = 2*ic - 1, j = jc, k = P.K_Min()-1; //k=0
          P_condensed(ic, jc, P_condensed.K_Min()-1) = P(i  ,j  ,k) 
            + P(i+1,j  ,k);
        }
    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int ic = P_condensed.I_Min(); ic <= P_condensed.I_Max(); ic++)
        for(int jc = P_condensed.J_Min(); jc <= P_condensed.J_Max(); jc++){
          int i = 2*ic - 1, j = jc, k = P.K_Max()+1; //k=kmax+1
          P_condensed(ic, jc, P_condensed.K_Max()+1) = P(i  ,j  ,k) 
            + P(i+1,j  ,k);
        } 
    }
}
//*****************************************************************************
// Extrapolates P_condensed to P (twice the size): coarse to fine
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Extend_Array(ARRAY_3D<T>& P_condensed, ARRAY_3D<T>& P)
{
  assert(P.Halo_Size()); assert(P_condensed.Halo_Size()); // halo > 0

  Fill_In_Eight_Boundary_Corners(P_condensed);

  if(!parameters->two_d){
    for(int ic = P_condensed.I_Min()-1; ic <= P_condensed.I_Max(); ic++)
      for(int jc = P_condensed.J_Min()-1; jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min()-1; kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic, j = 2*jc, k = 2*kc; 
          T p000=P_condensed(ic,jc,kc),     p100=P_condensed(ic+1,jc,kc),
            p010=P_condensed(ic,jc+1,kc),   p001=P_condensed(ic,jc,kc+1),
            p011=P_condensed(ic,jc+1,kc+1), p101=P_condensed(ic+1,jc,kc+1),
            p110=P_condensed(ic+1,jc+1,kc), p111=P_condensed(ic+1,jc+1,kc+1);
          // updating 8 corners of each sub-cell
          // Corner: 000
          P(i,j,k)   += (T).015625*((T)27*p000 + (T)9*(p100+p010+p001) 
              + (T)3*(p011+p101+p110) + p111);
          // Corner: 100
          P(i+1,j,k) += (T).015625*((T)27*p100 + (T)9*(p000+p110+p101) 
              + (T)3*(p111+p001+p010) + p011);
          // Corner: 010
          P(i,j+1,k) += (T).015625*((T)27*p010 + (T)9*(p110+p000+p011) 
              + (T)3*(p001+p111+p100) + p101);
          // Corner: 001
          P(i,j,k+1) += (T).015625*((T)27*p001 + (T)9*(p101+p011+p000) 
              + (T)3*(p010+p100+p111) + p110);
          // Corner: 011
          P(i,j+1,k+1) += (T).015625*((T)27*p011 + (T)9*(p111+p001+p010) 
              + (T)3*(p000+p110+p101) + p100 );
          // Corner: 101
          P(i+1,j,k+1) += (T).015625*((T)27*p101 + (T)9*(p001+p111+p100) 
              + (T)3*(p110+p000+p011) + p010);
          // Corner: 110
          P(i+1,j+1,k) += (T).015625*((T)27*p110 + (T)9*(p010+p100+p111) 
              + (T)3*(p101+p011+p000) + p001);
          // Corner: 111
          P(i+1,j+1,k+1) += (T).015625*((T)27*p111 + (T)9*(p011+p101+p110) 
              + (T)3*(p100+p010+p001)+p000);
        }//for_loops: ic,jc,kc
  } else {
    for(int ic = P_condensed.I_Min()-1; ic <= P_condensed.I_Max(); ic++)
      for(int jc = P_condensed.J_Min()-1; jc <= P_condensed.J_Max(); jc++)
        for(int kc = P_condensed.K_Min()-1; kc <= P_condensed.K_Max(); kc++){
          int i = 2*ic, j = jc, k = 2*kc; 
          /*
          T p000=P_condensed(ic,jc,kc),     p100=P_condensed(ic+1,jc,kc),
            p010=P_condensed(ic,jc+1,kc),   p001=P_condensed(ic,jc,kc+1),
            p011=P_condensed(ic,jc+1,kc+1), p101=P_condensed(ic+1,jc,kc+1),
            p110=P_condensed(ic+1,jc+1,kc), p111=P_condensed(ic+1,jc+1,kc+1);
          // updating 8 corners of each sub-cell
          // Corner: 000
          P(i,j,k)   += p000;
          // Corner: 100
          P(i+1,j,k) += p000; 
          // Corner: 001
          P(i,j,k+1) += p000;
          // Corner: 101
          P(i+1,j,k+1) += p000;
          */
          T p000=P_condensed(ic,jc,kc),   p100=P_condensed(ic+1,jc,kc),
            p001=P_condensed(ic,jc,kc+1), p101=P_condensed(ic+1,jc,kc+1);
          // updating 8 corners of each sub-cell
          // Corner: 000
          P(i,j,k)   += (T).015625*((T)27*p000 + (T)6*(p100+p001) 
              + (T)1*(p101) );
          // Corner: 100
          P(i+1,j,k) += (T).015625*((T)27*p100 + (T)6*(p000+p101) 
              + (T)1*(p001) );
          // Corner: 001
          P(i,j,k+1) += (T).015625*((T)27*p001 + (T)6*(p101+p000) 
              + (T)1*(p100) );
          // Corner: 101
          P(i+1,j,k+1) += (T).015625*((T)27*p101 + (T)6*(p001+p100) 
              + (T)1*(p000) );
          
        }//for_loops: ic,jc,kc
  }
}
//*****************************************************************************
// Fills in 8 boundary corners for P array
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Fill_In_Eight_Boundary_Corners(ARRAY_3D<T>& P)
{
  assert(P.Halo_Size()); // halo > 0

  int imin=P.I_Min(), imax=P.I_Max(), jmin=P.J_Min(), jmax=P.J_Max(), 
      kmin=P.K_Min(), kmax=P.K_Max();

  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
      if(mpi_driver->suth_proc == MPI_PROC_NULL)
        P(imin-1,jmin-1,kmin-1) = P(imin-1,jmin,kmin) + P(imin,jmin-1,kmin)
          + P(imin,jmin,kmin-1) - (T)2*P(imin,jmin,kmin);
      if(mpi_driver->nrth_proc == MPI_PROC_NULL)
        P(imin-1,jmin-1,kmax+1) = P(imin-1,jmin,kmax) + P(imin,jmin-1,kmax)
          + P(imin,jmin,kmax+1) - (T)2*P(imin,jmin,kmax);
    }//frnt

    if(mpi_driver->back_proc == MPI_PROC_NULL) {
      if(mpi_driver->suth_proc == MPI_PROC_NULL)
        P(imin-1,jmax+1,kmin-1) = P(imin-1,jmax,kmin) + P(imin,jmax+1,kmin)
          + P(imin,jmax,kmin-1) - (T)2*P(imin,jmax,kmin);
      if(mpi_driver->nrth_proc == MPI_PROC_NULL)
        P(imin-1,jmax+1,kmax+1) = P(imin-1,jmax,kmax) + P(imin,jmax+1,kmax)
          + P(imin,jmax,kmax+1) - (T)2*P(imin,jmax,kmax);
    }//back
  }//west

  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
      if(mpi_driver->suth_proc == MPI_PROC_NULL)
        P(imax+1,jmin-1,kmin-1) = P(imax+1,jmin,kmin) + P(imax,jmin-1,kmin)
          + P(imax,jmin,kmin-1) - (T)2*P(imax,jmin,kmin);
      if(mpi_driver->nrth_proc == MPI_PROC_NULL)
        P(imax+1,jmin-1,kmax+1) = P(imax+1,jmin,kmax) + P(imax,jmin-1,kmax)
          + P(imax,jmin,kmax+1) - (T)2*P(imax,jmin,kmax);
    }//frnt

    if(mpi_driver->back_proc == MPI_PROC_NULL) {
      if(mpi_driver->suth_proc == MPI_PROC_NULL)
        P(imax+1,jmax+1,kmin-1) = P(imax+1,jmax,kmin) + P(imax,jmax+1,kmin)
          + P(imax,jmax,kmin-1) - (T)2*P(imax,jmax,kmin);
      if(mpi_driver->nrth_proc == MPI_PROC_NULL)
        P(imax+1,jmax+1,kmax+1) = P(imax+1,jmax,kmax) + P(imax,jmax+1,kmax)
          + P(imax,jmax,kmax+1) - (T)2*P(imax,jmax,kmax);
    }//back
  }//east
}
//*****************************************************************************
// Main Smoothing function: used at every level of MG
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Smooth_Pressure(int level, 
    T& residual_l2, T& rhs_l2, T& residual_min, T& residual_max, 
    ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Residual)
{ 
  T residual_l2_old = (T)0;  
  METRIC_QUANTITIES<T> mq, mq_new; // metric quantities helper structures 
  // Initialize_P_RHS_Residual(level, P, RHS, Residual); // set based on level
  Initialize_Metric_Quantities(level, mq, mq_new); // set based on level
  
  // save initial residual in 'residual_l2_old'
  Compute_Residual(residual_l2_old, rhs_l2, residual_min, residual_max, 
      P, RHS, Residual, mq, mq_new);

  // exit if residual is zero
  if(!residual_l2_old) {residual_l2 = residual_l2_old; return;}

  // Main Smoothing Iteration Loop
  for(int iter=1; iter <= parameters->mg_max_smoothing_iters; iter++){

    // sub-iterations
    for(int n=1; n<=parameters->mg_smoothing_sub_iters; n++){
      //P(P.I_Min(),P.J_Min(),P.K_Min()) = 0.;
      Smooth_Pressure_In_Z(P, RHS, Residual, mq, mq_new); 

      //P(P.I_Min(),P.J_Min(),P.K_Min()) = 0.; 
      if(!parameters->two_d)
        Smooth_Pressure_In_Y(P, RHS, Residual, mq, mq_new); 

      //P(P.I_Min(),P.J_Min(),P.K_Min()) = 0.; 
      Smooth_Pressure_In_X(P, RHS, Residual, mq, mq_new); 

      //for(int i = P.I_Min()-1; i <= P.I_Max()+1; i++)
      //for(int k = P.K_Min()-1; k <= P.K_Max()+1; k++)
      // P(i,P.J_Min(),k) = 0.; 
      //P(P.I_Min(),P.J_Min(),P.K_Min()) = 0.; 
    }

    // check residual for convergence
    Compute_Residual(residual_l2, rhs_l2, residual_min, residual_max, 
        P, RHS, Residual, mq, mq_new);

    assert(residual_l2_old);
    if(residual_l2 / residual_l2_old <= parameters->mg_smoothing_converg_thresh)
      residual_l2_old = residual_l2;
    else break; // ends main iteration loop

  }//main loop
}
//*****************************************************************************
// Setup MQs based on the level: mq is a pointer structure, mq_new stores a copy
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Initialize_Metric_Quantities(int level, 
    METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new)
{
  // set pointers to metric quantities to be used as function parameters
  if(!level)
    mq.Set_Pointers(grid->G11, grid->G12, grid->G13, 
        grid->G21, grid->G22, grid->G23, 
        grid->G31, grid->G32, grid->G33, grid->GCC);
  else
    mq.Set_Pointers((*grid->G11_sub)(level), (*grid->G12_sub)(level), 
        (*grid->G13_sub)(level), (*grid->G21_sub)(level), 
        (*grid->G22_sub)(level), (*grid->G23_sub)(level), 
        (*grid->G31_sub)(level), (*grid->G32_sub)(level), 
        (*grid->G33_sub)(level), (*grid->GCC_sub)(level));

  // create copies of metric quantities (will be deleted with mq_new)
  mq_new.Create_Copies(mq.G11, mq.G12, mq.G13, mq.G21, mq.G22, mq.G23, 
      mq.G31, mq.G32, mq.G33, mq.GCC);
  mq_new.GCC->Set_All_Elements_To((T)0); // replacing GCC in mq_new

  int imin = mq_new.G11->I_Min(), imax = mq_new.G11->I_Max(),
      jmin = mq_new.G11->J_Min(), jmax = mq_new.G11->J_Max(),
      kmin = mq_new.G11->K_Min(), kmax = mq_new.G11->K_Max(), 
      halo = mq_new.G11->Halo_Size();  
  // boundary: left
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    int i = imin-1; //0
    for(int j = jmin-halo; j <= jmax+halo; j++)
      for(int k = kmin-halo; k <= kmax+halo; k++) {
        (*mq_new.G11)(i,j,k) = (T)0;
        (*mq_new.G12)(i,j,k) = (T)0;
        (*mq_new.G13)(i,j,k) = (T)0;
      }
  }
  // boundary: right
  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    int i = imax; //ii
    for(int j = jmin-halo; j <= jmax+halo; j++)
      for(int k = kmin-halo; k <= kmax+halo; k++) {
        (*mq_new.G11)(i,j,k) = (T)0;
        (*mq_new.G12)(i,j,k) = (T)0;
        (*mq_new.G13)(i,j,k) = (T)0;
      }
  }
  // boundary: front
  if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
    int j = jmin-1; //0
    for(int i = imin-halo; i <= imax+halo; i++)
      for(int k = kmin-halo; k <= kmax+halo; k++) {
        (*mq_new.G21)(i,j,k) = (T)0;
        (*mq_new.G22)(i,j,k) = (T)0;
        (*mq_new.G23)(i,j,k) = (T)0;
      }
  }
  // boundary: back
  if(mpi_driver->back_proc == MPI_PROC_NULL) {
    int j = jmax; //jj
    for(int i = imin-halo; i <= imax+halo; i++)
      for(int k = kmin-halo; k <= kmax+halo; k++) {
        (*mq_new.G21)(i,j,k) = (T)0;
        (*mq_new.G22)(i,j,k) = (T)0;
        (*mq_new.G23)(i,j,k) = (T)0;
      }
  }
  // boundary: south
  if(mpi_driver->suth_proc == MPI_PROC_NULL) {
    int k = kmin-1; //0
    for(int i = imin-halo; i <= imax+halo; i++)
      for(int j = jmin-halo; j <= jmax+halo; j++) {
        (*mq_new.G31)(i,j,k) = (T)0;
        (*mq_new.G32)(i,j,k) = (T)0;
        (*mq_new.G33)(i,j,k) = (T)0;
      }
  }
  // boundary: north
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
    int k = kmax; //kk
    for(int i = imin-halo; i <= imax+halo; i++)
      for(int j = jmin-halo; j <= jmax+halo; j++) {
        (*mq_new.G31)(i,j,k) = (T)0;
        (*mq_new.G32)(i,j,k) = (T)0;
        (*mq_new.G33)(i,j,k) = (T)0;
      }
  }
  // initialize new GCC
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++)
        (*mq_new.GCC)(i,j,k) = (*mq_new.G11)(i,j,k) + (*mq_new.G11)(i-1,j,k) 
          + (*mq_new.G22)(i,j,k) + (*mq_new.G22)(i,j-1,k) 
          + (*mq_new.G33)(i,j,k) + (*mq_new.G33)(i,j,k-1);
}
//*****************************************************************************
// Smoothing pressure in X Dimension
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Smooth_Pressure_In_X(ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, 
    ARRAY_3D<T>& Res, METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new)
{
  // LHS diagonals and RHS for tridiagonal linear system solve
  ARRAY_2D<T> A(P.J_Min(), P.J_Max(), P.I_Min()-1, P.I_Max()+1, 0),
    B(P.J_Min(), P.J_Max(), P.I_Min()-1, P.I_Max()+1, 0), 
    C(P.J_Min(), P.J_Max(), P.I_Min()-1, P.I_Max()+1, 0), 
    F(P.J_Min(), P.J_Max(), P.I_Min()-1, P.I_Max()+1, 0);    
  TRIDIAGONAL_SOLVER<T> LS_Solver(*mpi_driver);
  // Residual in the interior region
  for(int i = P.I_Min(); i <= P.I_Max(); i++)
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++) {
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq_new.G22)(i  ,j  ,k  ) * P(i  ,j+1,k  )
            + (*mq_new.G22)(i  ,j-1,k  ) * P(i  ,j-1,k  )
            + (*mq_new.G33)(i  ,j  ,k  ) * P(i  ,j  ,k+1)
            + (*mq_new.G33)(i  ,j  ,k-1) * P(i  ,j  ,k-1) );
        Res(i,j,k) += 
          (*mq_new.G12)(i  ,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i+1,j+1,k  ) - P(i+1,j-1,k  ))
          - (*mq_new.G12)(i-1,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i-1,j+1,k  ) - P(i-1,j-1,k  ))
          + (*mq_new.G13)(i  ,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i+1,j  ,k+1) - P(i+1,j  ,k-1))
          - (*mq_new.G13)(i-1,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i-1,j  ,k+1) - P(i-1,j  ,k-1));
        Res(i,j,k) +=   
          (*mq_new.G23)(i,j  ,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j+1,k+1) - P(i  ,j+1,k-1))
          - (*mq_new.G23)(i,j-1,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j-1,k+1) - P(i  ,j-1,k-1))
          + (*mq_new.G21)(i,j  ,k)*(P(i+1,j+1,k  ) - P(i-1,j+1,k  ))
          - (*mq_new.G21)(i,j-1,k)*(P(i+1,j-1,k  ) - P(i-1,j-1,k  ));
        Res(i,j,k) += 
          (*mq_new.G31)(i,j,k  )*(P(i+1,j  ,k+1) - P(i-1,j  ,k+1))
          - (*mq_new.G31)(i,j,k-1)*(P(i+1,j  ,k-1) - P(i-1,j  ,k-1))
          + (*mq_new.G32)(i,j,k  )*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k+1) - P(i  ,j-1,k+1))
          - (*mq_new.G32)(i,j,k-1)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k-1) - P(i  ,j-1,k-1));
      }

  Fill_In_Residual_On_Six_Boundary_Faces(Res, P, RHS, mq);

  // Setting up linear system
  for(int k = P.K_Min(); k <= P.K_Max(); k++) {

    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int j = P.J_Min(); j <= P.J_Max(); j++){
        A(j,i) = - (*mq_new.G11)(i-1,j,k)
          + (*mq_new.G21)(i,j,k) - (*mq_new.G21)(i,j-1,k) 
          + (*mq_new.G31)(i,j,k) - (*mq_new.G31)(i,j,k-1);
        C(j,i) = - (*mq_new.G11)(i  ,j,k) 
          - (*mq_new.G21)(i,j,k) + (*mq_new.G21)(i,j-1,k)
          - (*mq_new.G31)(i,j,k) + (*mq_new.G31)(i,j,k-1); 
        B(j,i) =   (*mq_new.GCC)(i,j,k);
        F(j,i) =   Res(i,j,k);
      }//for: i,j
    //bc: west
    if(mpi_driver->west_proc == MPI_PROC_NULL) {
      int i = P.I_Min()-1;
      for(int j = P.J_Min(); j <= P.J_Max(); j++){
        A(j,i) =   (T)0;
        B(j,i) =   (*mq.G11)(i,j,k);
        C(j,i) = - (*mq.G11)(i,j,k);
        F(j,i) =   Res(i,j,k);
      }
    }
    //bc: east
    if(mpi_driver->east_proc == MPI_PROC_NULL) {
      int i = P.I_Max()+1;
      for(int j = P.J_Min(); j <= P.J_Max(); j++){
        A(j,i) =   (*mq.G11)(i-1,j,k);
        B(j,i) = - (*mq.G11)(i-1,j,k);
        C(j,i) =   (T)0;
        F(j,i) =   Res(i,j,k);
      }
    }

    // solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(A, B, C, F, 
        parameters->periodic_in_x, mpi_driver->west_proc, mpi_driver->east_proc);

    // save solution F
    for(int i = Res.I_Min()-1; i <= Res.I_Max()+1; i++)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(i,j,k) = F(j,i);
  }//for: k

  // Fix other 4 boundaries
  //bc: frnt
  if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
    int j = Res.J_Min()-1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) = Res(i,j+1,k) + Res(i,j,k) / (*mq.G22)(i,j,k);
  }
  //bc: back
  if(mpi_driver->back_proc == MPI_PROC_NULL) {
    int j = Res.J_Max()+1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) =  Res(i,j-1,k) - Res(i,j,k) / (*mq.G22)(i,j-1,k);
  }
  //bc: south
  if(mpi_driver->suth_proc == MPI_PROC_NULL) {
    int k = Res.K_Min()-1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(i,j,k) =  Res(i,j,k+1) + Res(i,j,k) / (*mq.G33)(i,j,k);
  }
  //bc: north
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
    int k = Res.K_Max()+1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(i,j,k) =  Res(i,j,k-1) - Res(i,j,k) / (*mq.G33)(i,j,k-1);
  }

  Fill_In_Twelve_Boundary_Edges(Res);

  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(Res,1); 
  P = Res;
  //for(int i = P.I_Min()-1; i <= P.I_Max()+1; i++)
  // for(int j = P.J_Min()-1; j <= P.J_Max()+1; j++)
  //  for(int k = P.K_Min()-1; k <= P.K_Max()+1; k++)
  //    P(i,j,k) = Res(i,j,k);
}
//*****************************************************************************
// Smoothing pressure in Y Dimension
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Smooth_Pressure_In_Y(
    ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Res, 
    METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new)
{
  // LHS diagonals and RHS for tridiagonal linear system solve
  ARRAY_2D<T> A(P.I_Min(), P.I_Max(), P.J_Min()-1, P.J_Max()+1, 0),
    B(P.I_Min(), P.I_Max(), P.J_Min()-1, P.J_Max()+1, 0), 
    C(P.I_Min(), P.I_Max(), P.J_Min()-1, P.J_Max()+1, 0), 
    F(P.I_Min(), P.I_Max(), P.J_Min()-1, P.J_Max()+1, 0);    
  TRIDIAGONAL_SOLVER<T> LS_Solver(*mpi_driver);

  // Residual in the interior region
  for(int i = P.I_Min(); i <= P.I_Max(); i++)
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++) {
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq_new.G11)(i  ,j  ,k  ) * P(i+1,j  ,k  )
            + (*mq_new.G11)(i-1,j  ,k  ) * P(i-1,j  ,k  )
            + (*mq_new.G33)(i  ,j  ,k  ) * P(i  ,j  ,k+1)
            + (*mq_new.G33)(i  ,j  ,k-1) * P(i  ,j  ,k-1) );
        Res(i,j,k) += 
          (*mq_new.G12)(i  ,j,k)*(P(i+1,j+1,k  ) - P(i+1,j-1,k  ))
          - (*mq_new.G12)(i-1,j,k)*(P(i-1,j+1,k  ) - P(i-1,j-1,k  ))
          + (*mq_new.G13)(i  ,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i+1,j  ,k+1) - P(i+1,j  ,k-1))
          - (*mq_new.G13)(i-1,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i-1,j  ,k+1) - P(i-1,j  ,k-1));
        Res(i,j,k) +=   
          (*mq_new.G23)(i,j  ,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j+1,k+1) - P(i  ,j+1,k-1))
          - (*mq_new.G23)(i,j-1,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j-1,k+1) - P(i  ,j-1,k-1))
          + (*mq_new.G21)(i,j  ,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j+1,k  ) - P(i-1,j+1,k  ))
          - (*mq_new.G21)(i,j-1,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j-1,k  ) - P(i-1,j-1,k  ));
        Res(i,j,k) += 
          (*mq_new.G31)(i,j,k  )*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k+1) - P(i-1,j  ,k+1))
          - (*mq_new.G31)(i,j,k-1)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k-1) - P(i-1,j  ,k-1))
          + (*mq_new.G32)(i,j,k  )*(P(i  ,j+1,k+1) - P(i  ,j-1,k+1))
          - (*mq_new.G32)(i,j,k-1)*(P(i  ,j+1,k-1) - P(i  ,j-1,k-1));
      }

  Fill_In_Residual_On_Six_Boundary_Faces(Res, P, RHS, mq);

  // Setting up linear system
  for(int k = P.K_Min(); k <= P.K_Max(); k++) {

    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int j = P.J_Min(); j <= P.J_Max(); j++){
        A(i,j) = - (*mq_new.G22)(i,j-1,k)
          + (*mq_new.G12)(i,j,k) - (*mq_new.G12)(i-1,j,k) 
          + (*mq_new.G32)(i,j,k) - (*mq_new.G32)(i,j,k-1);
        C(i,j) = - (*mq_new.G22)(i  ,j,k) 
          - (*mq_new.G12)(i,j,k) + (*mq_new.G12)(i-1,j,k)
          - (*mq_new.G32)(i,j,k) + (*mq_new.G32)(i,j,k-1); 
        B(i,j) =   (*mq_new.GCC)(i,j,k);
        F(i,j) =   Res(i,j,k);
      }//for: i,j
    //bc: front
    if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
      int j = P.J_Min()-1;
      for(int i = P.I_Min(); i <= P.I_Max(); i++){
        A(i,j) =   (T)0;
        B(i,j) =   (*mq.G22)(i,j,k);
        C(i,j) = - (*mq.G22)(i,j,k);
        F(i,j) =   Res(i,j,k);
      }
    }
    //bc: back
    if(mpi_driver->back_proc == MPI_PROC_NULL) {
      int j = P.J_Max()+1;
      for(int i = P.I_Min(); i <= P.I_Max(); i++){
        A(i,j) =   (*mq.G22)(i,j-1,k);
        B(i,j) = - (*mq.G22)(i,j-1,k);
        C(i,j) =   (T)0;
        F(i,j) =   Res(i,j,k);
      }
    }

    // solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(A, B, C, F, 
        parameters->periodic_in_y, mpi_driver->frnt_proc, mpi_driver->back_proc);

    // save solution F
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int j = Res.J_Min()-1; j <= Res.J_Max()+1; j++)
        Res(i,j,k) = F(i,j);
  }//end_for: k

  // Fix the other 4 boundaries
  //bc: west
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    int i = Res.I_Min()-1;
    for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) = Res(i+1,j,k) + Res(i,j,k) / (*mq.G11)(i,j,k);
  }
  //bc: east
  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    int i = Res.I_Max()+1;
    for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) =  Res(i-1,j,k) - Res(i,j,k) / (*mq.G11)(i-1,j,k);
  }
  //bc: south
  if(mpi_driver->suth_proc == MPI_PROC_NULL) {
    int k = Res.K_Min()-1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(i,j,k) =  Res(i,j,k+1) + Res(i,j,k) / (*mq.G33)(i,j,k);
  }
  //bc: north
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
    int k = Res.K_Max()+1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(i,j,k) =  Res(i,j,k-1) - Res(i,j,k) / (*mq.G33)(i,j,k-1);
  }

  Fill_In_Twelve_Boundary_Edges(Res);

  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(Res,1); 
  P = Res;
  //for(int i = P.I_Min()-1; i <= P.I_Max()+1; i++)
  // for(int j = P.J_Min()-1; j <= P.J_Max()+1; j++)
  //  for(int k = P.K_Min()-1; k <= P.K_Max()+1; k++)      
  //    P(i,j,k) = Res(i,j,k);
}
//*****************************************************************************
// Smoothing pressure in Z Dimension
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Smooth_Pressure_In_Z(
    ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Res, 
    METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new)
{
  // LHS diagonals and RHS for tridiagonal linear system solve
  ARRAY_2D<T> A(P.I_Min(), P.I_Max(), P.K_Min()-1, P.K_Max()+1, 0),
    B(P.I_Min(), P.I_Max(), P.K_Min()-1, P.K_Max()+1, 0), 
    C(P.I_Min(), P.I_Max(), P.K_Min()-1, P.K_Max()+1, 0), 
    F(P.I_Min(), P.I_Max(), P.K_Min()-1, P.K_Max()+1, 0);    
  TRIDIAGONAL_SOLVER<T> LS_Solver(*mpi_driver);

  // Residual in the interior region
  for(int i = P.I_Min(); i <= P.I_Max(); i++)
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++) {
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq_new.G11)(i  ,j  ,k  ) * P(i+1,j  ,k  )
            + (*mq_new.G11)(i-1,j  ,k  ) * P(i-1,j  ,k  )
            + (*mq_new.G22)(i  ,j  ,k  ) * P(i  ,j+1,k  )
            + (*mq_new.G22)(i  ,j-1,k  ) * P(i  ,j-1,k  ) );
        Res(i,j,k) += 
          (*mq_new.G12)(i  ,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i+1,j+1,k  ) - P(i+1,j-1,k  ))
          - (*mq_new.G12)(i-1,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i-1,j+1,k  ) - P(i-1,j-1,k  ))
          + (*mq_new.G13)(i  ,j,k)*(P(i+1,j  ,k+1) - P(i+1,j  ,k-1))
          - (*mq_new.G13)(i-1,j,k)*(P(i-1,j  ,k+1) - P(i-1,j  ,k-1));
        Res(i,j,k) +=   
          (*mq_new.G23)(i,j  ,k)*(P(i  ,j+1,k+1) - P(i  ,j+1,k-1))
          - (*mq_new.G23)(i,j-1,k)*(P(i  ,j-1,k+1) - P(i  ,j-1,k-1))
          + (*mq_new.G21)(i,j  ,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j+1,k  ) - P(i-1,j+1,k  ))
          - (*mq_new.G21)(i,j-1,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j-1,k  ) - P(i-1,j-1,k  ));
        Res(i,j,k) += 
          (*mq_new.G31)(i,j,k  )*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k+1) - P(i-1,j  ,k+1))
          - (*mq_new.G31)(i,j,k-1)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k-1) - P(i-1,j  ,k-1))
          + (*mq_new.G32)(i,j,k  )*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k+1) - P(i  ,j-1,k+1))
          - (*mq_new.G32)(i,j,k-1)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k-1) - P(i  ,j-1,k-1));
      }

  Fill_In_Residual_On_Six_Boundary_Faces(Res, P, RHS, mq);

  // Setting up linear system
  for(int j = P.J_Min(); j <= P.J_Max(); j++){


    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++) {
        A(i,k) = - (*mq_new.G33)(i,j,k-1)
          + (*mq_new.G13)(i,j,k) - (*mq_new.G13)(i-1,j,k) 
          + (*mq_new.G23)(i,j,k) - (*mq_new.G23)(i,j-1,k);
        C(i,k) = - (*mq_new.G33)(i  ,j,k) 
          - (*mq_new.G13)(i,j,k) + (*mq_new.G13)(i-1,j,k)
          - (*mq_new.G23)(i,j,k) + (*mq_new.G23)(i,j-1,k); 
        B(i,k) =   (*mq_new.GCC)(i,j,k);
        F(i,k) =   Res(i,j,k);
      }//end_for: i,k

    //bc: south
    if(mpi_driver->suth_proc == MPI_PROC_NULL) {
      int k = P.K_Min()-1;
      for(int i = P.I_Min(); i <= P.I_Max(); i++){
        A(i,k) =   (T)0;
        B(i,k) =   (*mq.G33)(i,j,k);
        C(i,k) = - (*mq.G33)(i,j,k);
        F(i,k) =   Res(i,j,k);
      }
    }

    //bc: north
    if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
      int k = P.K_Max()+1;
      for(int i = P.I_Min(); i <= P.I_Max(); i++){
        A(i,k) =   (*mq.G33)(i,j,k-1);
        B(i,k) = - (*mq.G33)(i,j,k-1);
        C(i,k) =   (T)0;
        F(i,k) =   Res(i,j,k);
      }
    }

    // solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(A, B, C, F, 
        parameters->periodic_in_z, mpi_driver->suth_proc, mpi_driver->nrth_proc);

    // save solution F
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int k = Res.K_Min()-1; k <= Res.K_Max()+1; k++)
        Res(i,j,k) = F(i,k);
  }//end_for: k

  // Fix the other 4 boundaries
  //bc: west
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    int i = Res.I_Min()-1;
    for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) = Res(i+1,j,k) + Res(i,j,k) / (*mq.G11)(i,j,k);
  }

  //bc: east
  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    int i = Res.I_Max()+1;
    for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) =  Res(i-1,j,k) - Res(i,j,k) / (*mq.G11)(i-1,j,k);
  }

  //bc: front
  if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
    int j = Res.J_Min()-1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) =  Res(i,j+1,k) + Res(i,j,k) / (*mq.G22)(i,j,k);
  }

  //bc: back
  if(mpi_driver->back_proc == MPI_PROC_NULL) {
    int j = Res.J_Max()+1;
    for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(i,j,k) =  Res(i,j-1,k) - Res(i,j,k) / (*mq.G22)(i,j-1,k);
  }

  Fill_In_Twelve_Boundary_Edges(Res);

  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(Res,1); 
  P = Res; 
  //for(int i = P.I_Min()-1; i <= P.I_Max()+1; i++)
  // for(int j = P.J_Min()-1; j <= P.J_Max()+1; j++)
  //  for(int k = P.K_Min()-1; k <= P.K_Max()+1; k++)
  //    P(i,j,k) = Res(i,j,k);
}
//*****************************************************************************
// Adjusts residual values on 6 boundary edges. 
// Helper function: used by X,Y,Z smoothers.
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Fill_In_Residual_On_Six_Boundary_Faces(ARRAY_3D<T>& Res, 
    ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, METRIC_QUANTITIES<T>& mq)
{
  // Residual on left boundary
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    int i = P.I_Min()-1;
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq.G12)(i,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
                              + P(i+1,j+1,k  ) - P(i+1,j-1,k  ))
            + (*mq.G13)(i,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i+1,j  ,k+1) - P(i+1,j  ,k-1)) );
  }
  // Residual on right boundary
  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    int i = P.I_Max()+1;
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq.G12)(i-1,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
                                + P(i-1,j+1,k  ) - P(i-1,j-1,k  ))
            + (*mq.G13)(i-1,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i-1,j  ,k+1) - P(i-1,j  ,k-1)) );
  }
  // Residual on front boundary
  if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
    int j = P.J_Min()-1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq.G23)(i,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
                              + P(i  ,j+1,k+1) - P(i  ,j+1,k-1))
            + (*mq.G21)(i,j,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j+1,k  ) - P(i-1,j+1,k  )) );
  }
  // Residual on back boundary
  if(mpi_driver->back_proc == MPI_PROC_NULL) {
    int j = P.J_Max()+1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq.G23)(i,j-1,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
                                + P(i  ,j-1,k+1) - P(i  ,j-1,k-1))
            + (*mq.G21)(i,j-1,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j-1,k  ) - P(i-1,j-1,k  )) );
  }
  // Residual on south boundary
  if(mpi_driver->suth_proc == MPI_PROC_NULL) {
    int k = P.K_Min()-1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int j = P.J_Min(); j <= P.J_Max(); j++)
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq.G31)(i,j,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
                              + P(i+1,j  ,k+1) - P(i-1,j  ,k+1))
            + (*mq.G32)(i,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k+1) - P(i  ,j-1,k+1)) );
  }
  // Residual on north boundary
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
    int k = P.K_Max()+1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int j = P.J_Min(); j <= P.J_Max(); j++)
        Res(i,j,k) = - RHS(i,j,k) +
          ( (*mq.G31)(i,j,k-1)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
                                + P(i+1,j  ,k-1) - P(i-1,j  ,k-1))
            + (*mq.G32)(i,j,k-1)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k-1) - P(i  ,j-1,k-1)) );
  }
}
//*****************************************************************************
// Adjusts values of 12 outer edges of 3D array. 
// Helper function: used by X,Y,Z smoothers.
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Fill_In_Twelve_Boundary_Edges(ARRAY_3D<T>& Res)
{
  int imin = Res.I_Min(), imax = Res.I_Max(), jmin = Res.J_Min(), 
      jmax = Res.J_Max(), kmin = Res.K_Min(), kmax = Res.K_Max();
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(imin-1,jmin-1,k) = (T).5*(Res(imin,jmin-1,k) + Res(imin-1,jmin,k));

    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(imin-1,jmax+1,k) = (T).5*(Res(imin,jmax+1,k) + Res(imin-1,jmax,k));

    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(imin-1,j,kmin-1) = (T).5*(Res(imin,j,kmin-1) + Res(imin-1,j,kmin));

    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(imin-1,j,kmax+1) = (T).5*(Res(imin,j,kmax+1) + Res(imin-1,j,kmax));
  }//endif: west

  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(imax+1,jmin-1,k) = (T).5*(Res(imax,jmin-1,k) + Res(imax+1,jmin,k));

    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int k = Res.K_Min(); k <= Res.K_Max(); k++)
        Res(imax+1,jmax+1,k) = (T).5*(Res(imax,jmax+1,k) + Res(imax+1,jmax,k));

    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(imax+1,j,kmin-1) = (T).5*(Res(imax,j,kmin-1) + Res(imax+1,j,kmin));

    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int j = Res.J_Min(); j <= Res.J_Max(); j++)
        Res(imax+1,j,kmax+1) = (T).5*(Res(imax,j,kmax+1) + Res(imax+1,j,kmax));
  }//endif: east

  if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
        Res(i,jmin-1,kmin-1) = (T).5*(Res(i,jmin,kmin-1) + Res(i,jmin-1,kmin));

    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
        Res(i,jmin-1,kmax+1) = (T).5*(Res(i,jmin,kmax+1) + Res(i,jmin-1,kmax));
  }//endif: frnt

  if(mpi_driver->back_proc == MPI_PROC_NULL) {
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
        Res(i,jmax+1,kmin-1) = (T).5*(Res(i,jmax,kmin-1) + Res(i,jmax+1,kmin));

    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int i = Res.I_Min(); i <= Res.I_Max(); i++)
        Res(i,jmax+1,kmax+1) = (T).5*(Res(i,jmax,kmax+1) + Res(i,jmax+1,kmax));
  }//endif: back
}
//*****************************************************************************
// Computes residual of 'L(P) = RHS'
// Note: T**<-mq_new.G** and G**<-mq.G** in original implementation
//*****************************************************************************
  template<class T>
void PRESSURE<T>::Compute_Residual(T& residual_l2, T& rhs_l2, T& residual_min, 
    T& residual_max, ARRAY_3D<T>& P, ARRAY_3D<T>& RHS, ARRAY_3D<T>& Residual, 
    METRIC_QUANTITIES<T>& mq, METRIC_QUANTITIES<T>& mq_new)
{  
  assert(P.Halo_Size()); assert(RHS.Halo_Size());
  assert(Residual.Halo_Size()); // halo > 0

  T res_loc_l2=(T)0, rhs_loc_l2=(T)0, res_loc_min=(T)0, res_loc_max=(T)0;
  residual_l2 = rhs_l2 = residual_min = residual_max = (T)0;
  Residual.Set_All_Elements_To((T)0);

  // Residual in the interior region
  for(int i = P.I_Min(); i <= P.I_Max(); i++)
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++) {
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq_new.G11)(i  ,j  ,k  )*(P(i+1,j  ,k  ) - P(i,j,k))
            + (*mq_new.G11)(i-1,j  ,k  )*(P(i-1,j  ,k  ) - P(i,j,k))
            + (*mq_new.G22)(i  ,j  ,k  )*(P(i  ,j+1,k  ) - P(i,j,k))
            + (*mq_new.G22)(i  ,j-1,k  )*(P(i  ,j-1,k  ) - P(i,j,k))
            + (*mq_new.G33)(i  ,j  ,k  )*(P(i  ,j  ,k+1) - P(i,j,k))
            + (*mq_new.G33)(i  ,j  ,k-1)*(P(i  ,j  ,k-1) - P(i,j,k)) );
        Residual(i,j,k) -= 
          (*mq_new.G12)(i  ,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i+1,j+1,k  ) - P(i+1,j-1,k  ))
          - (*mq_new.G12)(i-1,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i-1,j+1,k  ) - P(i-1,j-1,k  ))
          + (*mq_new.G13)(i  ,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i+1,j  ,k+1) - P(i+1,j  ,k-1))
          - (*mq_new.G13)(i-1,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i-1,j  ,k+1) - P(i-1,j  ,k-1));
        Residual(i,j,k) -=   
          (*mq_new.G23)(i,j  ,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j+1,k+1) - P(i  ,j+1,k-1))
          - (*mq_new.G23)(i,j-1,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j-1,k+1) - P(i  ,j-1,k-1))
          + (*mq_new.G21)(i,j  ,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j+1,k  ) - P(i-1,j+1,k  ))
          - (*mq_new.G21)(i,j-1,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j-1,k  ) - P(i-1,j-1,k  ));
        Residual(i,j,k) -= 
          (*mq_new.G31)(i,j,k  )*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k+1) - P(i-1,j  ,k+1))
          - (*mq_new.G31)(i,j,k-1)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k-1) - P(i-1,j  ,k-1))
          + (*mq_new.G32)(i,j,k  )*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k+1) - P(i  ,j-1,k+1))
          - (*mq_new.G32)(i,j,k-1)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k-1) - P(i  ,j-1,k-1));
        res_loc_l2 += pow(Residual(i,j,k), 2);
        rhs_loc_l2 += pow(RHS(i,j,k), 2);
        res_loc_min = fmin(res_loc_min, Residual(i,j,k));
        res_loc_max = fmax(res_loc_max, Residual(i,j,k));
      }//endfor-loops:internal
  // Residual on left boundary
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    int i = P.I_Min()-1;
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq.G11)(i,j,k)*(P(i+1,j  ,k  ) - P(i  ,j  ,k  ))
            + (*mq.G12)(i,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i+1,j+1,k  ) - P(i+1,j-1,k  ))
            + (*mq.G13)(i,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i+1,j  ,k+1) - P(i+1,j  ,k-1)) );
  }
  // Residual on right boundary
  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    int i = P.I_Max()+1;
    for(int j = P.J_Min(); j <= P.J_Max(); j++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq.G11)(i-1,j,k)*(P(i  ,j  ,k  ) - P(i-1,j  ,k  ))
            + (*mq.G12)(i-1,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i-1,j+1,k  ) - P(i-1,j-1,k  ))
            + (*mq.G13)(i-1,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i-1,j  ,k+1) - P(i-1,j  ,k-1)) );
  }
  // Residual on front boundary
  if(mpi_driver->frnt_proc == MPI_PROC_NULL) {
    int j = P.J_Min()-1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq.G22)(i,j,k)*(P(i  ,j+1,k  ) - P(i  ,j  ,k  ))
            + (*mq.G23)(i,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j+1,k+1) - P(i  ,j+1,k-1))
            + (*mq.G21)(i,j,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j+1,k  ) - P(i-1,j+1,k  )) );
  }
  // Residual on back boundary
  if(mpi_driver->back_proc == MPI_PROC_NULL) {
    int j = P.J_Max()+1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int k = P.K_Min(); k <= P.K_Max(); k++)
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq.G22)(i,j-1,k)*(P(i  ,j  ,k  ) - P(i  ,j-1,k  ))
            + (*mq.G23)(i,j-1,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k-1)
              + P(i  ,j-1,k+1) - P(i  ,j-1,k-1))
            + (*mq.G21)(i,j-1,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j-1,k  ) - P(i-1,j-1,k  )) );
  }
  // Residual on south boundary
  if(mpi_driver->suth_proc == MPI_PROC_NULL) {
    int k = P.K_Min()-1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int j = P.J_Min(); j <= P.J_Max(); j++)
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq.G33)(i,j,k)*(P(i  ,j  ,k+1) - P(i  ,j  ,k  ))
            + (*mq.G31)(i,j,k)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k+1) - P(i-1,j  ,k+1))
            + (*mq.G32)(i,j,k)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k+1) - P(i  ,j-1,k+1)) );
  }
  // Residual on north boundary
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
    int k = P.K_Max()+1;
    for(int i = P.I_Min(); i <= P.I_Max(); i++)
      for(int j = P.J_Min(); j <= P.J_Max(); j++)
        Residual(i,j,k) = RHS(i,j,k) -
          ( (*mq.G33)(i,j,k-1)*(P(i  ,j  ,k  ) - P(i  ,j  ,k-1))
            + (*mq.G31)(i,j,k-1)*(P(i+1,j  ,k  ) - P(i-1,j  ,k  )
              + P(i+1,j  ,k-1) - P(i-1,j  ,k-1))
            + (*mq.G32)(i,j,k-1)*(P(i  ,j+1,k  ) - P(i  ,j-1,k  )
              + P(i  ,j+1,k-1) - P(i  ,j-1,k-1)) );
  }
  // collect results on root proc
  MPI_Comm& comm = mpi_driver->Get_Grid_Communicator();
  MPI_Reduce(&res_loc_l2,  &residual_l2,  1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&rhs_loc_l2,  &rhs_l2,       1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&res_loc_min, &residual_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&res_loc_max, &residual_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  // modify quantities on proc 0
  if(mpi_driver->my_rank == 0) {
      T denominator = P.I_Size() * mpi_driver->num_procs[0] * 
      P.J_Size() * mpi_driver->num_procs[1] * 
      P.K_Size() * mpi_driver->num_procs[2];
    assert(denominator);
    residual_l2 = sqrt(residual_l2 / denominator); 
    rhs_l2      = sqrt(rhs_l2 / denominator); 
  }
  // broacast results from root proc to everyone
  MPI_Bcast(&residual_l2,  1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&rhs_l2,       1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&residual_min, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&residual_max, 1, MPI_DOUBLE, 0, comm);
}
//*****************************************************************************
#endif
