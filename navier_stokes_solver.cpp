// navier_stokes_solver.cpp

#include "navier_stokes_solver.h"
#include <math.h>
//*****************************************************************************
// Constructor
//*****************************************************************************
  template<class T>
NAVIER_STOKES_SOLVER<T>::NAVIER_STOKES_SOLVER(int argc, char ** argv)
{
  //initialize simulation parameters class
  parameters = new PARAMETERS<T>(argc, argv);

  //initialize physical member variables: velocity, pressure, density(if needed)
  u = new ARRAY_3D<VECTOR_3D<T> >(parameters->i_min, parameters->i_max, 
      parameters->j_min, parameters->j_max, 
      parameters->k_min, parameters->k_max, parameters->halo_size);
  P = new ARRAY_3D<T>(parameters->i_min, parameters->i_max, 
      parameters->j_min, parameters->j_max, 
      parameters->k_min, parameters->k_max, parameters->halo_size-1);
  U_xi = new ARRAY_3D<T>(*P);
  U_et = new ARRAY_3D<T>(*P);
  U_zt = new ARRAY_3D<T>(*P);
  RHS_for_AB = new ARRAY_3D<VECTOR_3D<T> >(*P); //'hb' in Fcode

  //initialize scalar
  phi = new ARRAY_1D<ARRAY_3D<T>* >(parameters->num_scalars);
  if(parameters->scalar_advection){
    (*phi)(1) = new ARRAY_3D<T>(parameters->i_min, parameters->i_max, 
        parameters->j_min, parameters->j_max, 
        parameters->k_min, parameters->k_max, parameters->halo_size); //rho
    if(parameters->num_scalars==2)
      (*phi)(2) = new ARRAY_3D<T>(*(*phi)(1)); //passive scalar
  }
  else{ 
    (*phi)(1) = NULL;
  }

  //initialize mpi driver and grid (moving or stationary) classes
  mpi_driver = new MPI_DRIVER<T>(*parameters);
  if(parameters->moving_grid)
    grid = new CURVILINEAR_MOVING_GRID<T>(parameters, mpi_driver);
  else
    grid = new CURVILINEAR_GRID<T>(parameters, mpi_driver);
  turbulence = NULL; //turbulence is not implemented yet

  //initialize pressure, convection and scalar(E_p, if needed) classes
  pressure = new PRESSURE<T>(parameters, mpi_driver, grid, P, u,U_xi,U_et,U_zt);
  if(parameters->num_scalars <= 1){ //no scalar or density only
    convection1 = new CONVECTION<T>(parameters, mpi_driver, grid, (*phi)(1), 
        u, U_xi, U_et, U_zt); 
    convection2 = NULL;
  }
  else if (parameters->num_scalars == 2){ //density and passive scalar
    convection1 = new CONVECTION<T>(parameters, mpi_driver, grid, (*phi)(1), 
        u, U_xi, U_et, U_zt); 
    convection2 = new CONVECTION<T>(parameters, mpi_driver, grid, (*phi)(2), 
        u, U_xi, U_et, U_zt); 
  }
  if(parameters->scalar_advection){
    scalar1 = new SCALAR<T>(parameters,mpi_driver,grid,convection1,turbulence,(*phi)(1));
    if (parameters->num_scalars == 2) //density and passive scalar
      scalar2 = new SCALAR<T>(parameters,mpi_driver,grid,convection2,turbulence,(*phi)(2));
    else
      scalar2 = NULL;
    if(parameters->potential_energy) //depends on rho
      potential_energy = new POTENTIAL_ENERGY<T>(parameters,mpi_driver,grid,(*phi)(1));
    else 
      potential_energy = NULL;
  }
  else { 
    scalar1 = NULL;
    scalar2 = NULL;
  }

  //initialize moving grid engine (based on scalar)
  if(parameters->moving_grid && parameters->scalar_advection) 
    moving_grid_engine = new MOVING_GRID_ENGINE<T>(parameters, mpi_driver, 
        convection1, (CURVILINEAR_MOVING_GRID<T>*) grid, (*phi)(1), u);
  else moving_grid_engine = NULL;

  //aggregate any physical parameter in a timeseries
  if(parameters->aggregate_data) 
    data_aggregator = new DATA_AGGREGATOR<T>(parameters,mpi_driver,grid,u);

  //setup and save initial data as ts=0 or restart the sim
  if(!parameters->restart_timestep){
    //sets initial u&rho (requires MPI communication)
    Set_Initial_Conditions();
    //wait for all procs to get here and then save data for ts=0
    mpi_driver->Syncronize_All_Procs();
    //Save_Simulation_Data(); 
    if(mpi_driver->my_rank == 0) Save_Binary_Simulation_Parameters();
  }else Load_Simulation_Data_For_Restart(parameters->restart_timestep);
}
//*****************************************************************************
// Destructor
//*****************************************************************************
  template<class T>
NAVIER_STOKES_SOLVER<T>::~NAVIER_STOKES_SOLVER()
{
  delete u; delete U_xi; delete U_et; delete U_zt; delete P; delete RHS_for_AB;
  if(parameters->moving_grid) delete moving_grid_engine;
  if(parameters->scalar_advection){ 
    delete phi; 
    if(potential_energy) delete potential_energy;
    delete scalar1;
    if(parameters->num_scalars == 2) delete scalar2; 
  }
  if(parameters->aggregate_data) delete data_aggregator;
  delete grid; delete parameters; delete mpi_driver; 
}
//*****************************************************************************
// Predictor step: result is intermediate velocity (u*)
//*****************************************************************************
  template<class T>
void NAVIER_STOKES_SOLVER<T>::Predictor()
{
  T imin = grid->I_Min(), jmin = grid->J_Min(), kmin = grid->K_Min(), 
    imax = grid->I_Max(), jmax = grid->J_Max(), kmax = grid->K_Max(), 
    halo = grid->Halo_Size(), leading_term_AB_coefficient;
  ARRAY_3D<VECTOR_3D<T> > RHS(imin,imax,jmin,jmax,kmin,kmax, halo-1);//'su' in F
  TRIDIAGONAL_SOLVER<T> LS_Solver(*mpi_driver);
  ARRAY_2D<T> *A, *B, *C;       //LHS diagonals for tridiagonal solve
  ARRAY_2D<VECTOR_3D<T> > *F;   //RHS for tridiagonal solve

  /* Assemble the RHS */
  // Forward Euler timesteping on the first time step
  if(parameters->time_step == 1)
    leading_term_AB_coefficient = (T)1;
  else{
    leading_term_AB_coefficient = (T)1.5;
    // adding Adams-Bashforth contribution from previous time step
    RHS = (T)-.5 * (*RHS_for_AB);
  }

  // reset array for the new time step
  RHS_for_AB->Set_All_Elements_To(VECTOR_3D<T>((T)0,(T)0,(T)0));

  // add flow driving pressure gradient (+perturbation) term 
  if(parameters->pressure_gradient) Add_Pressure_Gradient_Term(RHS);

  // add convection term with QUICK
  convection1->Add_Quick_Scheme_Convection_Term(*RHS_for_AB);

  // add convection term to correct for moving grid velocity
  if(parameters->moving_grid)
    convection1->Add_Moving_Grid_Convection_Term(*RHS_for_AB);

  // LES self-similarity term
  if(parameters->turbulence) *RHS_for_AB -= *turbulence->tau;

  // Cross viscous terms at time = n-1 discretized with Crank-Nicolson
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++){
        //D_12 and D_13
        T visc_minus   = parameters->molecular_viscosity,
          visc_plus    = parameters->molecular_viscosity;
        VECTOR_3D<T>  visc_term;
        if(parameters->turbulence){
          visc_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i-1,j,k) );
          visc_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i+1,j,k) );
        }
        visc_term  = visc_plus * (*grid->G12)(i  ,j,k)*
          ((*u)(i,j+1,k) - (*u)(i,j-1,k) + (*u)(i+1,j+1,k) - (*u)(i+1,j-1,k));
        - visc_minus* (*grid->G12)(i-1,j,k)*
          ((*u)(i,j+1,k) - (*u)(i,j-1,k) + (*u)(i-1,j+1,k) - (*u)(i-1,j-1,k));
        visc_term += visc_plus * (*grid->G13)(i  ,j,k)*
          ((*u)(i,j,k+1) - (*u)(i,j,k-1) + (*u)(i+1,j,k+1) - (*u)(i+1,j,k-1));
        - visc_minus* (*grid->G13)(i-1,j,k)*
          ((*u)(i,j,k+1) - (*u)(i,j,k+1) + (*u)(i-1,j,k+1) - (*u)(i-1,j,k-1));

        //D_23 and D_21
        visc_minus = parameters->molecular_viscosity;
        visc_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          visc_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j-1,k) );
          visc_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j+1,k) );
        }
        visc_term += visc_plus * (*grid->G23)(i,j  ,k)*
          ((*u)(i,j,k+1) - (*u)(i,j,k-1) + (*u)(i,j+1,k+1) - (*u)(i,j+1,k-1));
        - visc_minus* (*grid->G23)(i,j-1,k)*
          ((*u)(i,j,k+1) - (*u)(i,j,k-1) + (*u)(i,j-1,k+1) - (*u)(i,j-1,k-1));
        visc_term += visc_plus * (*grid->G21)(i,j  ,k)*
          ((*u)(i+1,j,k) - (*u)(i-1,j,k) + (*u)(i+1,j+1,k) - (*u)(i-1,j+1,k));
        - visc_minus* (*grid->G21)(i,j-1,k)*
          ((*u)(i+1,j,k) - (*u)(i-1,j,k) + (*u)(i+1,j-1,k) - (*u)(i-1,j-1,k));

        //D_31 and D_32
        visc_minus = parameters->molecular_viscosity;
        visc_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          visc_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j,k-1) );
          visc_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j,k+1) );
        }
        visc_term += visc_plus * (*grid->G31)(i,j,k  )*
          ((*u)(i+1,j,k) - (*u)(i-1,j,k) + (*u)(i+1,j,k+1) - (*u)(i-1,j,k+1));
        - visc_minus* (*grid->G31)(i,j,k-1)*
          ((*u)(i+1,j,k) - (*u)(i-1,j,k) + (*u)(i+1,j,k-1) - (*u)(i-1,j,k-1));
        visc_term += visc_plus * (*grid->G32)(i,j,k  )*
          ((*u)(i,j+1,k) - (*u)(i,j-1,k) + (*u)(i,j+1,k+1) - (*u)(i,j-1,k+1));
        - visc_minus* (*grid->G32)(i,j,k-1)*
          ((*u)(i,j+1,k) - (*u)(i,j-1,k) + (*u)(i,j+1,k-1) - (*u)(i,j-1,k-1));

        // add D_12+D_13+D_23+D_21+D_31+D_32 to RHS
        (*RHS_for_AB)(i,j,k) += visc_term;
      }

  // add Coriolis and buoyancy source terms
  if(parameters->scalar_advection)
    for(int i = imin; i <= imax; i++)
      for(int j = jmin; j <= jmax; j++)
        for(int k = kmin; k <= kmax; k++){
          assert((*grid->inverse_Jacobian)(i,j,k));
          T jacobian = (T)1 / (*grid->inverse_Jacobian)(i,j,k);
          // Coriolis terms
          if(parameters->coriolis){
            (*RHS_for_AB)(i,j,k).x -=parameters->omega * (*u)(i,j,k).y*jacobian;
            (*RHS_for_AB)(i,j,k).y +=parameters->omega * (*u)(i,j,k).x*jacobian; 
          }
          // buoyancy term
          //(*RHS_for_AB)(i,j,k).z -= parameters->g * jacobian
          //                        * ((*rho)(i,j,k) - scalar->Rho_Rest(i,j,k));
          (*RHS_for_AB)(i,j,k).z -= parameters->g * jacobian
            * ((*(*phi)(1))(i,j,k) - scalar1->Rho_Rest(i,j,k));
        }

  // adding Adams-Bashforth contribution for current time step
  // 'RHS_for_AB' is saved from this point until the next time step
  //RHS += leading_term_AB_coefficient * (*RHS_for_AB);
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++)
        RHS(i,j,k) += leading_term_AB_coefficient * (*RHS_for_AB)(i,j,k);

  // Diagonal viscous terms at time = n-1 discretized with Crank-Nicolson
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++){
        //D_11
        T visc_minus   = parameters->molecular_viscosity,
          visc_plus    = parameters->molecular_viscosity;
        VECTOR_3D<T>  visc_term;
        if(parameters->turbulence){
          visc_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i-1,j,k) );
          visc_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i+1,j,k) );
        }
        visc_term=visc_plus * (*grid->G11)(i  ,j,k)*((*u)(i+1,j,k)-(*u)(i,j,k))
          - visc_minus* (*grid->G11)(i-1,j,k)*((*u)(i,j,k)-(*u)(i-1,j,k));
        //D_22
        visc_minus = parameters->molecular_viscosity;
        visc_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          visc_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j-1,k) );
          visc_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j+1,k) );
        }
        visc_term+=visc_plus* (*grid->G22)(i,j  ,k)*((*u)(i,j+1,k)-(*u)(i,j,k))
          - visc_minus* (*grid->G22)(i,j-1,k)*((*u)(i,j,k)-(*u)(i,j-1,k));
        //D_33
        visc_minus = parameters->molecular_viscosity;
        visc_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          visc_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j,k-1) );
          visc_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i,j,k+1) );
        }
        visc_term+=visc_plus* (*grid->G33)(i,j,k  )*((*u)(i,j,k+1)-(*u)(i,j,k))
          - visc_minus* (*grid->G33)(i,j,k-1)*((*u)(i,j,k)-(*u)(i,j,k-1));
        // add D_11+D_22+D_33 to RHS
        RHS(i,j,k) += visc_term;
      }

  // Jacobian difference adjustment term for a moving grid (component-wise prod)
  if(parameters->moving_grid)
    RHS -= (*((CURVILINEAR_MOVING_GRID<T>*)grid)->Jacobian_diff) * (*u);

  // Multiply every term by delta_time and inverse Jacobian
  // for(int i = RHS.I_Min_With_Halo(); i <= RHS.I_Max_With_Halo(); i++)
  //  for(int j = RHS.J_Min_With_Halo(); j <= RHS.J_Max_With_Halo(); j++)
  //    for(int k = RHS.K_Min_With_Halo(); k <= RHS.K_Max_With_Halo(); k++)
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++)
        RHS(i,j,k) *=  parameters->delta_time*(*grid->inverse_Jacobian)(i,j,k);

  /* END OF RHS ASSEMBLY */

  /* Solve 3 linear systems for momentum evolution eq.: M_i*M_j*M_k*x = RHS */
  //--------------------------------------------------------------------------
  // Solve for I-DIRECTION: (I-(dt/2J^{-1})*D_11)(u^*-u^n)_i = RHS
  //--------------------------------------------------------------------------
  A = new ARRAY_2D<T>(kmin,kmax, imin-1,imax+1, 0); //halo=0
  B = new ARRAY_2D<T>(kmin,kmax, imin-1,imax+1, 0);
  C = new ARRAY_2D<T>(kmin,kmax, imin-1,imax+1, 0);
  F = new ARRAY_2D<VECTOR_3D<T> >(kmin,kmax, imin-1,imax+1, 0);

  for(int j = jmin; j <= jmax; j++){
    // Construct tridiagonal matrices A,B,C for LHS and F for RHS(of system)
    for(int k = kmin; k <= kmax; k++)
      for(int i = imin; i <= imax; i++){
        T viscosity_minus = parameters->molecular_viscosity,
          viscosity_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          viscosity_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i-1,j,k) );
          viscosity_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k) + 
              (*turbulence->eddy_viscosity)(i+1,j,k) );
        }
        (*A)(k,i) = - viscosity_minus * (T).5 * parameters->delta_time * 
          (*grid->inverse_Jacobian)(i,j,k) * (*grid->G11)(i-1,j,k);
        (*C)(k,i) = - viscosity_plus  * (T).5 * parameters->delta_time * 
          (*grid->inverse_Jacobian)(i,j,k) * (*grid->G11)(i,j,k);
        (*B)(k,i) = (T)1 - (*A)(k,i) - (*C)(k,i);
        (*F)(k,i) = RHS(i,j,k);
      }// for: k,i
    // Boundary conditions for I-direction: west
    if(mpi_driver->west_proc == MPI_PROC_NULL)
      for(int k = kmin; k <= kmax; k++){
        (*A)(k,imin-1) = (T)0;	// imin-1 is 0 (in F code)
        (*B)(k,imin-1) = (T)1;	// free-slip: u_0 = u_1
        switch(parameters->west_bc){
          case FREE_SLIP: (*C)(k,imin-1) = (T)-1; break;
          case   NO_SLIP: (*C)(k,imin-1) = (T)1; break;
        }
        // RHS : BC west
        T dP_dXI = (*P)(imin  ,j  ,k  ) - (*P)(imin-1,j  ,k  ),
          dP_dET = (*P)(imin-1,j+1,k  ) - (*P)(imin-1,j-1,k  ) + 
            (*P)(imin  ,j+1,k  ) - (*P)(imin  ,j-1,k  ),
          dP_dZT = (*P)(imin-1,j  ,k+1) - (*P)(imin-1,j  ,k-1) + 
            (*P)(imin  ,j  ,k+1) - (*P)(imin  ,j  ,k-1);
        T coeff=(T)2*parameters->delta_time*(*grid->inverse_Jacobian)(imin,j,k);

        (*F)(k,imin-1).x = (*grid->XI_x)(imin-1,j,k) * dP_dXI + //imin-1||imin ?
          (T).125*((*grid->ET_x)(imin,j-1,k)+(*grid->ET_x)(imin,j,k)) * dP_dET +
          (T).125*((*grid->ZT_x)(imin,j,k-1)+(*grid->ZT_x)(imin,j,k)) * dP_dZT;
        (*F)(k,imin-1).y = (*grid->XI_y)(imin-1,j,k) * dP_dXI + //imin-1||imin ?
          (T).125*((*grid->ET_y)(imin,j-1,k)+(*grid->ET_y)(imin,j,k)) * dP_dET +
          (T).125*((*grid->ZT_y)(imin,j,k-1)+(*grid->ZT_y)(imin,j,k)) * dP_dZT;
        (*F)(k,imin-1).z = (*grid->XI_z)(imin-1,j,k) * dP_dXI + //imin-1||imin ?
          (T).125*((*grid->ET_z)(imin,j-1,k)+(*grid->ET_z)(imin,j,k)) * dP_dET +
          (T).125*((*grid->ZT_z)(imin,j,k-1)+(*grid->ZT_z)(imin,j,k)) * dP_dZT;
        (*F)(k,imin-1) *= coeff;
      }// west BC
    // Boundary conditions for I-direction: east
    if(mpi_driver->east_proc == MPI_PROC_NULL)
      for(int k = kmin; k <= kmax; k++){
        switch(parameters->east_bc){
          case FREE_SLIP: (*A)(k,imax+1) = (T)-1; break;// free-slip: u_0 = u_1
          case   NO_SLIP: (*A)(k,imax+1) = (T)1;  break;
        }
        (*B)(k,imax+1) = (T)1;
        (*C)(k,imax+1) = (T)0;
        // RHS : BC east
        T dP_dXI = (*P)(imax+1,j  ,k  ) - (*P)(imax  ,j  ,k  ),
          dP_dET = (*P)(imax  ,j+1,k  ) - (*P)(imax  ,j-1,k  ) + 
            (*P)(imax+1,j+1,k  ) - (*P)(imax+1,j-1,k  ),
          dP_dZT = (*P)(imax  ,j  ,k+1) - (*P)(imax  ,j  ,k-1) + 
            (*P)(imax+1,j  ,k+1) - (*P)(imax+1,j  ,k-1);
        T coeff=(T)2*parameters->delta_time*(*grid->inverse_Jacobian)(imax,j,k);

        (*F)(k,imax+1).x =(*grid-> XI_x)(imax,j,k) * dP_dXI + 
          (T).125*((*grid->ET_x)(imax,j-1,k)+(*grid->ET_x)(imax,j,k)) * dP_dET +
          (T).125*((*grid->ZT_x)(imax,j,k-1)+(*grid->ZT_x)(imax,j,k)) * dP_dZT;
        (*F)(k,imax+1).y = (*grid->XI_y)(imax,j,k) * dP_dXI + 
          (T).125*((*grid->ET_y)(imax,j-1,k)+(*grid->ET_y)(imax,j,k)) * dP_dET +
          (T).125*((*grid->ZT_y)(imax,j,k-1)+(*grid->ZT_y)(imax,j,k)) * dP_dZT;
        (*F)(k,imax+1).z = (*grid->XI_z)(imax,j,k) * dP_dXI + 
          (T).125*((*grid->ET_z)(imax,j-1,k)+(*grid->ET_z)(imax,j,k)) * dP_dET +
          (T).125*((*grid->ZT_z)(imax,j,k-1)+(*grid->ZT_z)(imax,j,k)) * dP_dZT;
        (*F)(k,imax+1) *= coeff;
      }//east BC

    //solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(*A, *B, *C, *F, 
        parameters->periodic_in_x, mpi_driver->west_proc, mpi_driver->east_proc);

    //solution of linear system is in F
    for(int k = kmin; k <= kmax; k++)
      for(int i = imin; i <= imax; i++)
        RHS(i,j,k) = (*F)(k,i);
  }//for: j
  delete A; delete B; delete C; delete F;

  //--------------------------------------------------------------------------
  // Solve for J-DIRECTION: (I-(dt/2J^{-1})*D_22)(u^*-u^n)_j = RHS 
  //--------------------------------------------------------------------------
  A = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0); //halo=0
  B = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0);
  C = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0);
  F = new ARRAY_2D<VECTOR_3D<T> >(imin,imax, jmin-1,jmax+1, 0);

  for(int k = kmin; k <= kmax; k++){
    // Construct tridiagonal matrices A,B,C for LHS and F for RHS (of LinSys)
    for(int i = imin; i <= imax; i++)
      for(int j = jmin; j <= jmax; j++){
        T viscosity_minus = parameters->molecular_viscosity,
          viscosity_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          viscosity_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j  ,k) + 
              (*turbulence->eddy_viscosity)(i,j-1,k) );
          viscosity_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j  ,k) + 
              (*turbulence->eddy_viscosity)(i,j+1,k) );
        }
        (*A)(i,j) = - viscosity_minus * (T).5 * parameters->delta_time * 
          (*grid->inverse_Jacobian)(i,j,k) * (*grid->G22)(i,j-1,k);
        (*C)(i,j) = - viscosity_plus  * (T).5 * parameters->delta_time * 
          (*grid->inverse_Jacobian)(i,j,k) * (*grid->G22)(i,j,k);
        (*B)(i,j) = (T)1 - (*A)(i,j) - (*C)(i,j);
        (*F)(i,j) = RHS(i,j,k);
      }// for: i,j
    // Boundary conditions for J-direction: front
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
        (*A)(i,jmin-1) = (T)0;	// jmin-1 is 0 (in F code)
        (*B)(i,jmin-1) = (T)1;	// no-slip: u_0 = -u_1
        switch(parameters->frnt_bc){
          case FREE_SLIP: (*C)(i,jmin-1) = (T)-1; break;
          case   NO_SLIP: (*C)(i,jmin-1) = (T)1;  break; // no-slip: u_0 = -u_1
        }
        // RHS : front BC
        T dP_dET = (*P)(i  ,jmin  ,k) - (*P)(i  ,jmin-1,k),
          dP_dXI = (*P)(i+1,jmin-1,k) - (*P)(i-1,jmin-1,k) + 
            (*P)(i+1,jmin  ,k) - (*P)(i-1,jmin  ,k),
          dP_dZT = (*P)(i,jmin-1,k+1) - (*P)(i,jmin-1,k-1) + 
            (*P)(i,jmin  ,k+1) - (*P)(i,jmin  ,k-1);
        T coeff=(T)2*parameters->delta_time*(*grid->inverse_Jacobian)(i,jmin,k);

        (*F)(i,jmin-1).x = 
          (T).125*((*grid->XI_x)(i-1,jmin,k)+(*grid->XI_x)(i,jmin,k)) * dP_dXI +
          /*jmin-1 or jmin?*/              (*grid->ET_x)(i,jmin-1,k)  * dP_dET +
          (T).125*((*grid->ZT_x)(i,jmin,k-1)+(*grid->ZT_x)(i,jmin,k)) * dP_dZT;
        (*F)(i,jmin-1).y = 
          (T).125*((*grid->XI_y)(i-1,jmin,k)+(*grid->XI_y)(i,jmin,k)) * dP_dXI +
          /*jmin-1 or jmin?*/              (*grid->ET_y)(i,jmin-1,k)  * dP_dET +
          (T).125*((*grid->ZT_y)(i,jmin,k-1)+(*grid->ZT_y)(i,jmin,k)) * dP_dZT;
        (*F)(i,jmin-1).z = 
          (T).125*((*grid->XI_z)(i-1,jmin,k)+(*grid->XI_z)(i,jmin,k)) * dP_dXI +
          /*jmin-1 or jmin?*/              (*grid->ET_z)(i,jmin-1,k)  * dP_dET +
          (T).125*((*grid->ZT_z)(i,jmin,k-1)+(*grid->ZT_z)(i,jmin,k)) * dP_dZT;
        (*F)(i,jmin-1) *= coeff;
      }// front BC
    // Boundary conditions for J-direction: back
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
        switch(parameters->back_bc){
          case FREE_SLIP: (*A)(i,jmax+1) = (T)-1; break; //free-slip:u_N = u_N+1
          case   NO_SLIP: (*A)(i,jmax+1) = (T)1;  break; 
        }	
        (*B)(i,jmax+1) = (T)1;
        (*C)(i,jmax+1) = (T)0;
        // RHS : BC
        T dP_dET = (*P)(i  ,jmax+1,k  ) - (*P)(i  ,jmax  ,k  ),
          dP_dXI = (*P)(i+1,jmax  ,k  ) - (*P)(i-1,jmax  ,k  ) + 
            (*P)(i+1,jmax+1,k  ) - (*P)(i-1,jmax+1,k  ),
          dP_dZT = (*P)(i  ,jmax  ,k+1) - (*P)(i  ,jmax  ,k-1) + 
            (*P)(i  ,jmax+1,k+1) - (*P)(i  ,jmax+1,k-1);
        T coeff=(T)2*parameters->delta_time*(*grid->inverse_Jacobian)(i,jmax,k);

        (*F)(i,jmax+1).x =
          (T).125*((*grid->XI_x)(i-1,jmax,k)+(*grid->XI_x)(i,jmax,k)) * dP_dXI +
          (*grid->ET_x)(i,jmax,k)  * dP_dET +
          (T).125*((*grid->ZT_x)(i,jmax,k-1)+(*grid->ZT_x)(i,jmax,k)) * dP_dZT;
        (*F)(i,jmax+1).y =
          (T).125*((*grid->XI_y)(i-1,jmax,k)+(*grid->XI_y)(i,jmax,k)) * dP_dXI +
          (*grid->ET_y)(i,jmax,k)  * dP_dET +
          (T).125*((*grid->ZT_y)(i,jmax,k-1)+(*grid->ZT_y)(i,jmax,k)) * dP_dZT;
        (*F)(i,jmax+1).z =
          (T).125*((*grid->XI_z)(i-1,jmax,k)+(*grid->XI_z)(i,jmax,k)) * dP_dXI +
          (*grid->ET_z)(i,jmax,k)  * dP_dET +
          (T).125*((*grid->ZT_z)(i,jmax,k-1)+(*grid->ZT_z)(i,jmax,k)) * dP_dZT;
        (*F)(i,jmax+1) *= coeff;
      }//back BC

    //solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(*A, *B, *C, *F, 
        parameters->periodic_in_y, mpi_driver->frnt_proc, mpi_driver->back_proc);

    //solution of linear system is in F
    for(int i = imin; i <= imax; i++)
      for(int j = jmin; j <= jmax; j++)
        RHS(i,j,k) = (*F)(i,j);
  }//for: k
  delete A; delete B; delete C; delete F;  

  //--------------------------------------------------------------------------
  // Solve for K-DIRECTION: (I-(dt/2J^{-1})*D_33)(u^*-u^n)_k = RHS
  //--------------------------------------------------------------------------
  A = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0); //halo=0
  B = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0);
  C = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0);
  F = new ARRAY_2D<VECTOR_3D<T> >(imin,imax, kmin-1,kmax+1, 0);

  for(int j = jmin; j <= jmax; j++){
    // Construct tridiagonal matrices A,B,C for LHS and F for RHS (of LinSys)
    for(int k = kmin; k <= kmax; k++)
      for(int i = imin; i <= imax; i++){
        T viscosity_minus = parameters->molecular_viscosity,
          viscosity_plus  = parameters->molecular_viscosity;
        if(parameters->turbulence){
          viscosity_minus += (T).5*( (*turbulence->eddy_viscosity)(i,j,k  ) + 
              (*turbulence->eddy_viscosity)(i,j,k-1) );
          viscosity_plus  += (T).5*( (*turbulence->eddy_viscosity)(i,j,k  ) + 
              (*turbulence->eddy_viscosity)(i,j,k+1) );
        }
        (*A)(i,k) = - viscosity_minus * (T).5 * parameters->delta_time * 
          (*grid->inverse_Jacobian)(i,j,k) * (*grid->G33)(i,j,k-1);
        (*C)(i,k) = - viscosity_plus  * (T).5 * parameters->delta_time * 
          (*grid->inverse_Jacobian)(i,j,k) * (*grid->G33)(i,j,k);
        (*B)(i,k) = (T)1 - (*A)(i,k) - (*C)(i,k);
        (*F)(i,k) = RHS(i,j,k);
      }// for: i,k
    // Boundary conditions for K-direction: bottom
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
        (*A)(i,kmin-1) = (T)0;	// kmin-1 is 0 (in F code)
        (*B)(i,kmin-1) = (T)1;
        switch(parameters->suth_bc){
          case FREE_SLIP: (*C)(i,kmin-1) = (T)-1; break;// free-slip: u_0 = u_1
          case   NO_SLIP: (*C)(i,kmin-1) = (T)1;  break; 
        }
        // RHS : BC bottom
        T dP_dZT = (*P)(i  ,j  ,kmin  ) - (*P)(i  ,j  ,kmin-1),
          dP_dXI = (*P)(i+1,j  ,kmin-1) - (*P)(i-1,j  ,kmin-1) + 
            (*P)(i+1,j  ,kmin  ) - (*P)(i-1,j  ,kmin  ),
          dP_dET = (*P)(i  ,j+1,kmin-1) - (*P)(i  ,j-1,kmin-1) + 
            (*P)(i  ,j+1,kmin  ) - (*P)(i  ,j-1,kmin  );
        T coeff=(T)2*parameters->delta_time*(*grid->inverse_Jacobian)(i,j,kmin);

        (*F)(i,kmin-1).x = 
          (T).125*((*grid->XI_x)(i-1,j,kmin)+(*grid->XI_x)(i,j,kmin)) * dP_dXI +
          (T).125*((*grid->ET_x)(i,j-1,kmin)+(*grid->ET_x)(i,j,kmin)) * dP_dET +
          /*kmin-1 or kmin ? */             (*grid->ZT_x)(i,j,kmin-1) * dP_dZT;
        (*F)(i,kmin-1).y = 
          (T).125*((*grid->XI_y)(i-1,j,kmin)+(*grid->XI_y)(i,j,kmin)) * dP_dXI +
          (T).125*((*grid->ET_y)(i,j-1,kmin)+(*grid->ET_y)(i,j,kmin)) * dP_dET +
          /*kmin-1 or kmin ? */             (*grid->ZT_y)(i,j,kmin-1) * dP_dZT;
        (*F)(i,kmin-1).z = 
          (T).125*((*grid->XI_z)(i-1,j,kmin)+(*grid->XI_z)(i,j,kmin)) * dP_dXI +
          (T).125*((*grid->ET_z)(i,j-1,kmin)+(*grid->ET_z)(i,j,kmin)) * dP_dET +
          /*kmin-1 or kmin ? */            (*grid->ZT_z)(i,j,kmin-1)  * dP_dZT;
        (*F)(i,kmin-1) *= coeff;
      }// bottom BC
    // Boundary conditions for K-direction: top
    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
        switch(parameters->nrth_bc){
          case FREE_SLIP: (*A)(i,kmax+1) = (T)-1; break;//free-slip: u_N = u_N+1
          case   NO_SLIP: (*A)(i,kmax+1) = (T)1;  break; 
        }
        (*B)(i,kmax+1) = (T)1;
        (*C)(i,kmax+1) = (T)0;
        // RHS : BC top
        T dP_dZT = (*P)(i  ,j  ,kmax+1) - (*P)(i  ,j  ,kmax  ),
          dP_dXI = (*P)(i+1,j  ,kmax  ) - (*P)(i-1,j  ,kmax  ) + 
            (*P)(i+1,j  ,kmax+1) - (*P)(i-1,j  ,kmax+1),
          dP_dET = (*P)(i  ,j+1,kmax  ) - (*P)(i  ,j-1,kmax  ) + 
            (*P)(i  ,j+1,kmax+1) - (*P)(i  ,j-1,kmax+1);
        T coeff=(T)2*parameters->delta_time*(*grid->inverse_Jacobian)(i,j,kmax);

        (*F)(i,kmax+1).x = 
          (T).125*((*grid->XI_x)(i-1,j,kmax)+(*grid->XI_x)(i,j,kmax)) * dP_dXI +
          (T).125*((*grid->ET_x)(i,j-1,kmax)+(*grid->ET_x)(i,j,kmax)) * dP_dET +
          (*grid->ZT_x)(i,j,kmax)  * dP_dZT;
        (*F)(i,kmax+1).y = 
          (T).125*((*grid->XI_y)(i-1,j,kmax)+(*grid->XI_y)(i,j,kmax)) * dP_dXI +
          (T).125*((*grid->ET_y)(i,j-1,kmax)+(*grid->ET_y)(i,j,kmax)) * dP_dET +
          (*grid->ZT_y)(i,j,kmax)  * dP_dZT;
        (*F)(i,kmax+1).z = 
          (T).125*((*grid->XI_z)(i-1,j,kmax)+(*grid->XI_z)(i,j,kmax)) * dP_dXI +
          (T).125*((*grid->ET_z)(i,j-1,kmax)+(*grid->ET_z)(i,j,kmax)) * dP_dET +
          (*grid->ZT_z)(i,j,kmax)  * dP_dZT;
        (*F)(i,kmax+1) *= coeff;
      }//top BC
    //solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(*A, *B, *C, *F, 
        parameters->periodic_in_z, mpi_driver->suth_proc, mpi_driver->nrth_proc);
    //solution of linear system is in F
    for(int i = imin; i <= imax; i++)
      for(int k = kmin; k <= kmax; k++)
        RHS(i,j,k) = (*F)(i,k);
  }//for: k
  delete A; delete B; delete C; delete F;

  //mpi_driver->Write_Global_Array_To_Disk("rhs", RHS,parameters->time_step);

  // Update U to be intermediate velocity (u*)
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++)
        (*u)(i,j,k) += RHS(i,j,k);

  // BCs for u*
  Linear_Extrapolate_Into_Halo_Regions(*u);
  //Enforce_Velocity_BC(*u);

  // Swap halo regions among procs
  mpi_driver->Exchange_Ghost_Values_For_Vector_Field(*u);  

  convection1->Quick_Velocity_Flux_Update(*u); // Velocity fluxes on faces

  mpi_driver->Syncronize_All_Procs(); // Wait for all procs to finish

  /*
  //test print
  T Qp = 0.;
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    for(int j = jmax; j >= jmin; j--) {
      cout << (*U_xi)(mpi_driver->local_grid_lower_bound[0]-1,j,ceil(kmax/2)) << endl;
      Qp += (*U_xi)(mpi_driver->local_grid_lower_bound[0]-1,j,ceil(kmax/2)); 
    }
    cout << "Qp = " << Qp << endl;
  }
  */
}
//*****************************************************************************
// Intermediate velocity correction due to enforce incompressibility
//*****************************************************************************
  template<class T>
void NAVIER_STOKES_SOLVER<T>::Corrector()
{
  // Correct velocity
  for(int i=grid->I_Min(); i<=grid->I_Max(); i++)
    for(int j=grid->J_Min(); j<=grid->J_Max(); j++)
      for(int k=grid->K_Min(); k<=grid->K_Max(); k++) {
        T P_east = (*P)(i+1,j,k) + (*P)(i,j,k), 
          P_west = (*P)(i-1,j,k) + (*P)(i,j,k),
          P_back = (*P)(i,j+1,k) + (*P)(i,j,k),
          P_frnt = (*P)(i,j-1,k) + (*P)(i,j,k),
          P_nrth = (*P)(i,j,k+1) + (*P)(i,j,k),
          P_suth = (*P)(i,j,k-1) + (*P)(i,j,k),
          coeff = (T).5*parameters->delta_time*(*grid->inverse_Jacobian)(i,j,k);

        (*u)(i,j,k).x -= coeff * 
          ( (*grid->XI_x)(i,j,k)*P_east - (*grid->XI_x)(i-1,j,k)*P_west
            + (*grid->ET_x)(i,j,k)*P_back - (*grid->ET_x)(i,j-1,k)*P_frnt
            + (*grid->ZT_x)(i,j,k)*P_nrth - (*grid->ZT_x)(i,j,k-1)*P_suth );

        (*u)(i,j,k).y -= coeff * 
          ( (*grid->XI_y)(i,j,k)*P_east - (*grid->XI_y)(i-1,j,k)*P_west
            + (*grid->ET_y)(i,j,k)*P_back - (*grid->ET_y)(i,j-1,k)*P_frnt
            + (*grid->ZT_y)(i,j,k)*P_nrth - (*grid->ZT_y)(i,j,k-1)*P_suth );

        (*u)(i,j,k).z -= coeff * 
          ( (*grid->XI_z)(i,j,k)*P_east - (*grid->XI_z)(i-1,j,k)*P_west
            + (*grid->ET_z)(i,j,k)*P_back - (*grid->ET_z)(i,j-1,k)*P_frnt
            + (*grid->ZT_z)(i,j,k)*P_nrth - (*grid->ZT_z)(i,j,k-1)*P_suth );

      }

  // Boundary Conditions and communication between processors
  Set_Progressive_Wave_BC(parameters->time);
  Enforce_Velocity_BC(*u);
  mpi_driver->Exchange_Ghost_Values_For_Vector_Field(*u);
  int *lower_boundary = mpi_driver->local_grid_lower_bound,
      *upper_boundary = mpi_driver->local_grid_upper_bound;

  // Correct velocity fluxes on faces
  for(int i=lower_boundary[0]; i<=upper_boundary[0]; i++)
    for(int j=grid->J_Min(); j<=grid->J_Max(); j++)
      for(int k=grid->K_Min(); k<=grid->K_Max(); k++) 
        (*U_xi)(i,j,k) -= parameters->delta_time *
          ( (*grid->G11)(i,j,k)*((*P)(i+1,j  ,k  ) - (*P)(i  ,j  ,k  ))
            +(*grid->G12)(i,j,k)*((*P)(i  ,j+1,k  ) - (*P)(i  ,j-1,k  )
              +(*P)(i+1,j+1,k  ) - (*P)(i+1,j-1,k  ))
            +(*grid->G13)(i,j,k)*((*P)(i  ,j  ,k+1) - (*P)(i  ,j  ,k-1)
              +(*P)(i+1,j  ,k+1) - (*P)(i+1,j  ,k-1)) );

  for(int i=grid->I_Min(); i<=grid->I_Max(); i++)
    for(int j=lower_boundary[1]; j<=upper_boundary[1]; j++)
      for(int k=grid->K_Min(); k<=grid->K_Max(); k++) 
        (*U_et)(i,j,k) -= parameters->delta_time *
          ( (*grid->G22)(i,j,k)*((*P)(i  ,j+1,k  ) - (*P)(i  ,j  ,k  ))
            +(*grid->G23)(i,j,k)*((*P)(i  ,j  ,k+1) - (*P)(i  ,j  ,k-1)
              +(*P)(i  ,j+1,k+1) - (*P)(i  ,j+1,k-1))
            +(*grid->G21)(i,j,k)*((*P)(i+1,j  ,k  ) - (*P)(i-1,j  ,k  )
              +(*P)(i+1,j+1,k  ) - (*P)(i-1,j+1,k  )) );

  for(int i=grid->I_Min(); i<=grid->I_Max(); i++)
    for(int j=grid->J_Min(); j<=grid->J_Max(); j++)
      for(int k=lower_boundary[2]; k<=upper_boundary[2]; k++)
        (*U_zt)(i,j,k) -= parameters->delta_time *
          ( (*grid->G33)(i,j,k)*((*P)(i  ,j  ,k+1) - (*P)(i  ,j  ,k  ))
            +(*grid->G31)(i,j,k)*((*P)(i+1,j  ,k  ) - (*P)(i-1,j  ,k  )
              +(*P)(i+1,j  ,k+1) - (*P)(i-1,j  ,k+1))
            +(*grid->G32)(i,j,k)*((*P)(i  ,j+1,k  ) - (*P)(i  ,j-1,k  )
              +(*P)(i  ,j+1,k+1) - (*P)(i  ,j-1,k+1)) );

  /*
  //test print
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    for(int j = grid->J_Max(); j >= grid->J_Min(); j--) {
      //cout << ( (*u)(grid->I_Min(),j,ceil(grid->K_Max()/2)).x + 
      //          (*u)(grid->I_Min()-1,j,ceil(grid->K_Max()/2)).x ) / 2. << endl;
      cout << (*u)(grid->I_Min()-2,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()-1,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()  ,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()+1,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()+2,j,ceil(grid->K_Max()/2)).x << "\t" <<
              ( (*u)(grid->I_Min(),j,ceil(grid->K_Max()/2)).x + 
                (*u)(grid->I_Min()-1,j,ceil(grid->K_Max()/2)).x ) / 2. << endl;
    }
  }
  */

  // Wait for all procs to finish
  mpi_driver->Syncronize_All_Procs();
}
//*****************************************************************************
// Extrapolating (linear with 3pts) vector field into the Halo regions
// Extrapolation rule for 3 directions: (0)<-(1,2,3), (-1)<-(0,1,2)
//*****************************************************************************
  template<class T>
void NAVIER_STOKES_SOLVER<T>::Linear_Extrapolate_Into_Halo_Regions(
    ARRAY_3D<VECTOR_3D<T> >& u)
{
  assert(u.Halo_Size() == 2);
  int xmin = u.I_Min(), xmax = u.I_Max(), ymin = u.J_Min(), ymax = u.J_Max(),
      zmin = u.K_Min(), zmax = u.K_Max(), halo = u.Halo_Size();

  if(mpi_driver->west_proc == MPI_PROC_NULL)
    for(int j = ymin-halo; j <= ymax+halo; j++)
      for(int k = zmin-halo; k <= zmax+halo; k++){
        u(xmin-1,j,k) = (T)3*(u(xmin,j,k)-u(xmin+1,j,k)) + u(xmin+2,j,k);
        u(xmin-2,j,k) = (T)3*(u(xmin-1,j,k)-u(xmin,j,k)) + u(xmin+1,j,k);   
      }

  if(mpi_driver->east_proc == MPI_PROC_NULL)
    for(int j = ymin-halo; j <= ymax+halo; j++)
      for(int k = zmin-halo; k <= zmax+halo; k++){
        u(xmax+1,j,k) = (T)3*(u(xmax,j,k)-u(xmax-1,j,k)) + u(xmax-2,j,k);
        u(xmax+2,j,k) = (T)3*(u(xmax+1,j,k)-u(xmax,j,k)) + u(xmax-1,j,k);  
      }

  if(mpi_driver->frnt_proc == MPI_PROC_NULL)
    for(int i = xmin-halo; i <= xmax+halo; i++)
      for(int k = zmin-halo; k <= zmax+halo; k++){
        u(i,ymin-1,k) = (T)3*(u(i,ymin,k)-u(i,ymin+1,k)) + u(i,ymin+2,k);
        u(i,ymin-2,k) = (T)3*(u(i,ymin-1,k)-u(i,ymin,k)) + u(i,ymin+1,k);   
      }

  if(mpi_driver->back_proc == MPI_PROC_NULL)
    for(int i = xmin-halo; i <= xmax+halo; i++)
      for(int k = zmin-halo; k <= zmax+halo; k++){
        u(i,ymax+1,k) = (T)3*(u(i,ymax,k)-u(i,ymax-1,k)) + u(i,ymax-2,k);
        u(i,ymax+2,k) = (T)3*(u(i,ymax+1,k)-u(i,ymax,k)) + u(i,ymax-1,k); 
      }

  if(mpi_driver->suth_proc == MPI_PROC_NULL)
    for(int i = xmin-halo; i <= xmax+halo; i++)
      for(int j = ymin-halo; j <= ymax+halo; j++){
        u(i,j,zmin-1) = (T)3*(u(i,j,zmin)-u(i,j,zmin+1)) + u(i,j,zmin+2);
        u(i,j,zmin-2) = (T)3*(u(i,j,zmin-1)-u(i,j,zmin)) + u(i,j,zmin+1);  
      }

  if(mpi_driver->nrth_proc == MPI_PROC_NULL)
    for(int i = xmin-halo; i <= xmax+halo; i++)
      for(int j = ymin-halo; j <= ymax+halo; j++){
        u(i,j,zmax+1) = (T)3*(u(i,j,zmax)-u(i,j,zmax-1)) + u(i,j,zmax-2);
        u(i,j,zmax+2) = (T)3*(u(i,j,zmax+1)-u(i,j,zmax)) + u(i,j,zmax-1);     
      }
}
//*****************************************************************************
// Global velocity BC
//*****************************************************************************
  template<class T>
void NAVIER_STOKES_SOLVER<T>::Enforce_Velocity_BC(ARRAY_3D<VECTOR_3D<T> >& u)
{
  // BC: Top
  if(mpi_driver->nrth_proc == MPI_PROC_NULL)
    for(int i=u.I_Min_With_Halo(); i<=u.I_Max_With_Halo(); i++)
      for(int j=u.J_Min_With_Halo(); j<=u.J_Max_With_Halo(); j++) {
        if(parameters->lid_velocity){
          //u(i,j,u.K_Max()) = (*parameters->lid_velocity)(i,j);
          u(i,j,u.K_Max()+1) = (T)2 * (*parameters->lid_velocity)(i,j) - 
            u(i,j,u.K_Max());
        }
        else
          switch(parameters->nrth_bc){
            case FREE_SLIP: u(i,j,u.K_Max()+1) =   u(i,j,u.K_Max()); break;
            case   NO_SLIP: u(i,j,u.K_Max()+1) = - u(i,j,u.K_Max()); break;
          }
        u(i,j,u.K_Max()+2) = (T)3 * (u(i,j,u.K_Max()+1) - u(i,j,u.K_Max()))
          +  u(i,j,u.K_Max()-1);
      }
  // BC: Bottom
  if(mpi_driver->suth_proc == MPI_PROC_NULL)
    for(int i=u.I_Min_With_Halo(); i<=u.I_Max_With_Halo(); i++)
      for(int j=u.J_Min_With_Halo(); j<=u.J_Max_With_Halo(); j++) {
        if(parameters->bed_velocity){
          //u(i,j,u.K_Min()) = (*parameters->bed_velocity)(i,j);
          u(i,j,u.K_Min()-1) = (T)2 * (*parameters->bed_velocity)(i,j) - 
            u(i,j,u.K_Min());
        }
        else	  
          switch(parameters->suth_bc){
            case FREE_SLIP: u(i,j,u.K_Min()-1) =   u(i,j,u.K_Min()); break;
            case   NO_SLIP: u(i,j,u.K_Min()-1) = - u(i,j,u.K_Min()); break; 
          }
        u(i,j,u.K_Min()-2) = (T)3 * (u(i,j,u.K_Min()-1) - u(i,j,u.K_Min()))
          +  u(i,j,u.K_Min()+1);

        if(parameters->bed_velocity) 
          u(i,j,u.K_Min()-2) =  u(i,j,u.K_Min()-1) =
            (*parameters->bed_velocity)(i,j);
      }
  // BC: Left
  if(mpi_driver->west_proc == MPI_PROC_NULL)
    for(int j=u.J_Min_With_Halo(); j<=u.J_Max_With_Halo(); j++)
      for(int k=u.K_Min_With_Halo(); k<=u.K_Max_With_Halo(); k++) {
        if(parameters->west_velocity){
          u(u.I_Min()-1,j,k).x = (T)2 * (*parameters->west_velocity)(j,k).x - 
            u(u.I_Min(),j,k).x;                      //Set velocity in x
          u(u.I_Min()-1,j,k).y = u(u.I_Min(),j,k).y; //Gradient free in y,z
          u(u.I_Min()-1,j,k).z = u(u.I_Min(),j,k).z;
        }
        else
          switch(parameters->west_bc){
            case FREE_SLIP: u(u.I_Min()-1,j,k) =   u(u.I_Min(),j,k); break;
            case   NO_SLIP: u(u.I_Min()-1,j,k) = - u(u.I_Min(),j,k); break; 
          }
        u(u.I_Min()-2,j,k) = (T)3 * (u(u.I_Min()-1,j,k) - u(u.I_Min(),j,k))
          +  u(u.I_Min()+1,j,k);
      }
  // BC: Right
  if(mpi_driver->east_proc == MPI_PROC_NULL)
    for(int j=u.J_Min_With_Halo(); j<=u.J_Max_With_Halo(); j++)
      for(int k=u.K_Min_With_Halo(); k<=u.K_Max_With_Halo(); k++) {
        switch(parameters->east_bc){
          case FREE_SLIP: u(u.I_Max()+1,j,k) =   u(u.I_Max(),j,k); break;
          case   NO_SLIP: u(u.I_Max()+1,j,k) = - u(u.I_Max(),j,k); break; 
        }
        u(u.I_Max()+2,j,k) = (T)3 * (u(u.I_Max()+1,j,k) - u(u.I_Max(),j,k))
          +  u(u.I_Max()-1,j,k);
      }
  // BC: Back
  if(mpi_driver->back_proc == MPI_PROC_NULL)
    for(int i=u.I_Min_With_Halo(); i<=u.I_Max_With_Halo(); i++) 
      for(int k=u.K_Min_With_Halo(); k<=u.K_Max_With_Halo(); k++) {
        switch(parameters->back_bc){
          case FREE_SLIP: u(i,u.J_Max()+1,k) =   u(i,u.J_Max(),k); break;
          case   NO_SLIP: u(i,u.J_Max()+1,k) = - u(i,u.J_Max(),k); break; 
        }
        u(i,u.J_Max()+2,k) = (T)3 * (u(i,u.J_Max()+1,k) - u(i,u.J_Max(),k))
          +  u(i,u.J_Max()-1,k);

      }
  // BC: Front
  if(mpi_driver->frnt_proc == MPI_PROC_NULL)
    for(int i=u.I_Min_With_Halo(); i<=u.I_Max_With_Halo(); i++) 
      for(int k=u.K_Min_With_Halo(); k<=u.K_Max_With_Halo(); k++) {
        switch(parameters->frnt_bc){
          case FREE_SLIP: u(i,u.J_Min()-1,k) =   u(i,u.J_Min(),k); break;
          case   NO_SLIP: u(i,u.J_Min()-1,k) = - u(i,u.J_Min(),k); break; 
        }
        u(i,u.J_Min()-2,k) = (T)3 * (u(i,u.J_Min()-1,k) - u(i,u.J_Min(),k))
          +  u(i,u.J_Min()+1,k);
      }
}
//*****************************************************************************
// Calculates of CFL for the current time step and updates 'max_cfl'
// Note: value of CFL cannot get above 'critical_cfl'
//*****************************************************************************
  template<class T> 
T NAVIER_STOKES_SOLVER<T>::Calculate_CFL()
{
  T cfl = (T)0;
  // calculate local CFL
  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
        T cell_cfl;
        if(parameters->moving_grid)
          cell_cfl = 
            fabs( (*U_xi)(i-1,j,k)
                - (*((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_xi)(i-1,j,k)
                + (*U_xi)(i,j,k)
                - (*((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_xi)(i,j,k) )
            + fabs( (*U_et)(i,j-1,k)
                - (*((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_et)(i,j-1,k)
                + (*U_et)(i,j,k)
                - (*((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_et)(i,j,k) )
            + fabs( (*U_zt)(i,j,k-1) 
                - (*((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_zt)(i,j,k-1)
                + (*U_zt)(i,j,k)
                - (*((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_zt)(i,j,k) );
        else
          cell_cfl = fabs((*U_xi)(i-1,j,k) + (*U_xi)(i,j,k)) +
            fabs((*U_et)(i,j-1,k) + (*U_et)(i,j,k)) +
            fabs((*U_zt)(i,j,k-1) + (*U_zt)(i,j,k));
        cell_cfl *= (*grid->inverse_Jacobian)(i,j,k);
        cfl = fmax(cfl, fabs(cell_cfl));
      }
  cfl *= (T).5 * parameters->delta_time; // .5 is taken outside of for-loops

  // find MAX CFL among all procs and replace local CFL with MAX
  mpi_driver->Replace_With_Max_Value_Among_All_Procs(cfl);

  //keep track of the global MAX CFL
  parameters->max_cfl = fmax(parameters->max_cfl,cfl); 
  return cfl;
}
//*****************************************************************************
// Adds a forcing term associated with the Pressue Gradient (+ perturbation)
//*****************************************************************************
  template<class T> 
void NAVIER_STOKES_SOLVER<T>::Add_Pressure_Gradient_Term(
    ARRAY_3D<VECTOR_3D<T> >& RHS)
{
  T omega = parameters->freq_p_grad,
    P0 = parameters->amp_p_grad*omega*omega;
  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
        assert((*grid->inverse_Jacobian)(i,j,k));
        RHS(i,j,k)   += (*parameters->pressure_gradient);
        //RHS(i,j,k).z += P0*cos(omega*parameters->time); //perturbation
        RHS(i,j,k)   /= (*grid->inverse_Jacobian)(i,j,k);
      }
}
//*****************************************************************************
// Saving data aggregaged on Proc #0
//*****************************************************************************
  template<class T>
void NAVIER_STOKES_SOLVER<T>::Save_Simulation_Data()
{
  if(parameters->time_step==0)
    mpi_driver->Write_Global_Array_To_Disk("node_position", *grid->grid, 0, true);
  if(parameters->save_instant_velocity)
    mpi_driver->Write_Global_Array_To_Disk("velocity",*u,parameters->time_step, true);
  if(parameters->save_pressure)
    mpi_driver->Write_Global_Array_To_Disk("pressure",*P,parameters->time_step);
  if(parameters->scalar_advection){
    mpi_driver->Write_Global_Array_To_Disk("density", *(*phi)(1), 
        parameters->time_step);
    if(parameters->num_scalars == 2)
      mpi_driver->Write_Global_Array_To_Disk("scalar", *(*phi)(2), 
          parameters->time_step);
    if(parameters->potential_energy) potential_energy->Write_To_Disk();
  }  
  if(parameters->aggregate_data) data_aggregator->Write_To_Disk();
  if(parameters->save_fluxes){
    mpi_driver->Write_Global_Array_To_Disk("u_xi", *U_xi,parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("u_et", *U_et,parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("u_zt", *U_zt,parameters->time_step);
  }
  if(parameters->moving_grid){
    if(parameters->time_step) //for ts=0, already saved
      mpi_driver->Write_Global_Array_To_Disk("node_position", *grid->grid,
          parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("u_grid_xi", 
        *((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_xi,parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("u_grid_et", 
        *((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_et,parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("u_grid_zt", 
        *((CURVILINEAR_MOVING_GRID<T>*)grid)->U_grid_zt,parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("jac_diff", 
        *((CURVILINEAR_MOVING_GRID<T>*)grid)->Jacobian_diff,parameters->time_step);
    mpi_driver->Write_Global_Array_To_Disk("inverse_jacobian", 
        *grid->inverse_Jacobian,parameters->time_step);
    ((CURVILINEAR_MOVING_GRID<T>*)grid)->Write_Volume_Evolution_To_Disk();
  }
  //if(parameters->time_step==0){
  //  mpi_driver->Write_Binary_Local_Array("bin_p",*P);  
  //  mpi_driver->Read_Binary_Local_Array("bin_p",*P); 
  //}
 }
//*****************************************************************************
// Saving binary simulation data to output directory for each proc
//*****************************************************************************
  template<class T>
int NAVIER_STOKES_SOLVER<T>::Save_Binary_Simulation_Data()
{
  if(parameters->time_step==parameters->save_data_timestep_period)
    Save_Binary_Simulation_Data("grid", *grid->grid);
  if(parameters->save_instant_velocity)
    Save_Binary_Simulation_Data("velocity", *u);
  if(parameters->scalar_advection){
    Save_Binary_Simulation_Data("density", *(*phi)(1));
    if(parameters->num_scalars==2)
      Save_Binary_Simulation_Data("scalar", *(*phi)(2));
  }
  if(parameters->save_pressure)
    Save_Binary_Simulation_Data("pressure", *P);

  return 1;
}
//*****************************************************************************
// Helper: Saving binary output data to output directory for each proc
//*****************************************************************************
  template<class T>
int NAVIER_STOKES_SOLVER<T>::Save_Binary_Simulation_Data(string a_name, 
    ARRAY_3D<T>& a)
{
  stringstream filename;

  //overwrite existing file
  if(parameters->time_step==parameters->save_data_timestep_period){
    filename << parameters->output_dir << a_name << "."<< mpi_driver->my_rank;
    ofstream output(filename.str().c_str(), ios::out | ios::trunc | ios::binary);
    if(!output){
      cout << "ERROR: could not open "<< a_name << " file for writing" <<endl;
      return 0;
    }
    mpi_driver->Write_Binary_Local_Array_Output(output, a); 
    output.close();
  } 

  //append at end of existing file
  else{
    filename << parameters->output_dir << a_name << "."<< mpi_driver->my_rank;
    ofstream output(filename.str().c_str(), ios::out | ios::app | ios::binary);
    if(!output){
      cout << "ERROR: could not open "<< a_name << " file for writing" <<endl;
      return 0;
    }
    mpi_driver->Write_Binary_Local_Array_Output(output, a); 
    output.close();
  }

  return 1;
}
//*****************************************************************************
// Helper: Saving binary output data to output directory for each proc
//*****************************************************************************
  template<class T>
int NAVIER_STOKES_SOLVER<T>::Save_Binary_Simulation_Data(string a_name, 
    ARRAY_3D<VECTOR_3D<T> >& a)
{
  stringstream filename;

  //overwrite existing file
  if(parameters->time_step==parameters->save_data_timestep_period){
    filename << parameters->output_dir << a_name << "."<< mpi_driver->my_rank;
    ofstream output(filename.str().c_str(), ios::out | ios::trunc | ios::binary);
    if(!output){
      cout << "ERROR: could not open "<< a_name << " file for writing" <<endl;
      return 0;
    }
    mpi_driver->Write_Binary_Local_Array_Output(output, a); 
    output.close();
  } 

  //append at end of existing file
  else{
    filename << parameters->output_dir << a_name << "."<< mpi_driver->my_rank;
    ofstream output(filename.str().c_str(), ios::out | ios::app | ios::binary);
    if(!output){
      cout << "ERROR: could not open "<< a_name << " file for writing" <<endl;
      return 0;
    }
    mpi_driver->Write_Binary_Local_Array_Output(output, a); 
    output.close();
  }

  return 1;
}
//*****************************************************************************
// Output simulation parameters in binary
//*****************************************************************************
template<class T>
int NAVIER_STOKES_SOLVER<T>::Save_Binary_Simulation_Parameters()
{
  if(!mpi_driver->my_rank){
    stringstream filename;

    filename << parameters->output_dir << "parameters.0";
    ofstream output(filename.str().c_str(), ios::out | ios::trunc | ios::binary);
    if(!output){
      cout<<"ERROR: could not open parameters file for writing"<<endl;
      return 0;
    }

    output.write(reinterpret_cast<char *>(&parameters->num_total_nodes_x),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->num_total_nodes_y),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->num_total_nodes_z),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->num_cpu_x),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->num_cpu_y),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->num_cpu_z),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->max_timestep),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->save_data_timestep_period),sizeof(int));
    output.write(reinterpret_cast<char *>(&parameters->delta_time),sizeof(T));
    output.write(reinterpret_cast<char *>(&parameters->x_length),sizeof(T));
    output.write(reinterpret_cast<char *>(&parameters->y_length),sizeof(T));
    output.write(reinterpret_cast<char *>(&parameters->z_length),sizeof(T));
    output.close();
  }
  
  return 1;
}
//*****************************************************************************
// Saving binary restart data to output directory for each proc
//*****************************************************************************
  template<class T>
int NAVIER_STOKES_SOLVER<T>::Save_Simulation_Data_For_Restart()
{
  stringstream filename;
  filename << parameters->output_dir << "restart_t" 
    << parameters->time_step<< "."<< mpi_driver->my_rank;

  ofstream output(filename.str().c_str(), ios::out | ios::binary);
  if(!output){
    cout<<"ERROR: could not open RESTART file for writing"<<endl;
    return 0;
  }
  //parameters: time, time_step
  output.write(reinterpret_cast<char *>(&parameters->time),sizeof(T));
  output.write(reinterpret_cast<char *>(&parameters->time_step),sizeof(T));

  //TODO: grid: only for moving grid

  //physical parameters: u (U_xi/et/zt), P, rho(if convection)
  mpi_driver->Write_Binary_Local_Array(output, *u); 
  mpi_driver->Write_Binary_Local_Array(output, *U_xi); 
  mpi_driver->Write_Binary_Local_Array(output, *U_et); 
  mpi_driver->Write_Binary_Local_Array(output, *U_zt); 
  mpi_driver->Write_Binary_Local_Array(output, *P); 
  mpi_driver->Write_Binary_Local_Array(output, *RHS_for_AB); 
  output.close();
  return 1;
}
//*****************************************************************************
// Loading binary restart data from ./restart directory for each proc
//*****************************************************************************
  template<class T>
int NAVIER_STOKES_SOLVER<T>::Load_Simulation_Data_For_Restart(int restart_ts)
{
  stringstream filename;
  filename << parameters->output_dir << "restart_t" 
    << restart_ts<< "."<< mpi_driver->my_rank;
  ifstream input(filename.str().c_str(), ios::in | ios::binary);
  if(!input){
    cout<<"ERROR: could not open RESTART file for reading"<<endl;
    return 0;
  }
  //parameters: time, timestep
  input.read(reinterpret_cast<char *>(&parameters->time),sizeof(T));
  input.read(reinterpret_cast<char *>(&parameters->time_step),sizeof(T));
  //TODO: grid: only for moving grid
  //physical parameters: u (U_xi/et/zt), P, rho(if convection)
  mpi_driver->Read_Binary_Local_Array(input, *u);
  mpi_driver->Read_Binary_Local_Array(input, *U_xi);
  mpi_driver->Read_Binary_Local_Array(input, *U_et); 
  mpi_driver->Read_Binary_Local_Array(input, *U_zt);
  mpi_driver->Read_Binary_Local_Array(input, *P);
  mpi_driver->Read_Binary_Local_Array(input, *RHS_for_AB);
  input.close();
  return 1;
}
//*****************************************************************************
// Set initial conditions for physical variables (u and rho)
//*****************************************************************************
  template<class T> 
void NAVIER_STOKES_SOLVER<T>::Set_Initial_Conditions()
{
  //set Velocity profile
  //KH billows example
  /*
     for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++)
     for(int j=grid->J_Min_With_Halo(); j<=grid->J_Max_With_Halo(); j++)
     for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++){	
     T lambda = (T).25 * parameters->x_length,
     k_lambda = (T)2 * parameters->pi / lambda,
     zeta = (T).001 * sin( k_lambda * (*grid)(i,j,k).x ),
     temp = tanh( ((*grid)(i,j,k).y - parameters->y_length/2. + zeta)
     /(4.*parameters->y_length/parameters->num_total_nodes_y) );
     (*u)(i,j,k) = VECTOR_3D<T>(temp,(T)0,(T)0);
     (*rho)(i,j,k) = temp;
     }
     convection->Quick_Velocity_Flux_Update(*u); // velocity fluxes on faces
     */

  /*
  //set initial Density profile
  //Lock-exchange example
  for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++)
  for(int j=grid->J_Min_With_Halo(); j<=grid->J_Max_With_Halo(); j++)
  for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++)
  if( abs((*grid)(i,j,k).x) < (T).5*parameters->x_length ) 
  (*rho)(i,j,k) = (T).0005;
  else (*rho)(i,j,k) = (T)-.0005;
  */

  /*
  //Sloshing wave example
  //if(parameters->scalar_advection) scalar->Set_Initial_Density_Profile();
  if(parameters->scalar_advection)
  for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++)
  for(int j=grid->J_Min_With_Halo(); j<=grid->J_Max_With_Halo(); j++)
  for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++){
  if( abs((*grid)(i,j,k).x) > (T).45*parameters->x_length && 
  abs((*grid)(i,j,k).x) < (T).55*parameters->x_length &&
  abs((*grid)(i,j,k).y) > (T).45*parameters->y_length && 
  abs((*grid)(i,j,k).y) < (T).55*parameters->y_length) 
  (*rho)(i,j,k) = (T).0005;
  else (*rho)(i,j,k) = (T)-.0005;
  if( i==grid->I_Max() && 
  abs((*grid)(i,j,k).y) > (T).65*parameters->y_length && 
  abs((*grid)(i,j,k).y) < (T).75*parameters->y_length)
  (*rho)(i,j,k) = (T).0005;
  //(*u)(i,j,k) = VECTOR_3D<T>(1.,(T)0,(T)0);
  }
  */

  /*
  //Channel flow
  //srand(time(NULL)*(mpi_driver->my_rank+1)); //to prevent repeating patterns
  T kappa = .41;
  T u_star = .005;
  for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++)
  for(int j=grid->J_Min(); j<=grid->J_Max(); j++)
  for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++){
//(*u)(i,j,k).x = .1 * sqrt(1.+(*grid)(i,j,k).y/parameters->y_length);
T d = parameters->depth ? (*parameters->depth)(i,k) 
: parameters->y_length;
//  if((*grid)(i,j,k).y+d <= (T)0)
//  (*u)(i,j,k).x = 0.;
//else
(*u)(i,j,k).x = u_star/kappa * log(9. * u_star
   * ((*grid)(i,j,k).y+d)
   / parameters->molecular_viscosity);
// (*u)(i,j,k)  += .05 * (*u)(i,j,k).x * rand() / (T)RAND_MAX;
}
*/

//Variable depth channel example
//if(parameters->scalar_advection) scalar->Set_Uniform_Density_Profile((T)1);

//Bobby's cases
//Set velocity and stratification
if(parameters->scalar_advection) {
  T alpha = parameters->alpha;
  T delta = parameters->delta;
  T ratio = parameters->ratio; //delta_rho/rho_0
  T rho0 = parameters->rho0;;
  T delta_rho = ratio*rho0;  
  T interface_loc = parameters->z_length/parameters->upper_layer_depth;
  T delta_perturb = parameters->delta_perturb;
  T lambda_perturb = parameters->lambda_perturb;
  T a = parameters->a;
  T Lw = parameters->Lw;
  T zeta;

  if(parameters->solitary_wave){
    //Single Solitary wave 
    //density
    for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++) {
      for(int j=grid->J_Min_With_Halo(); j<=grid->J_Max_With_Halo(); j++) {
        for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++) {
          zeta = -a*exp(-pow((*grid)(i,j,k).x/Lw,2))
                  + delta_perturb*cos(2*parameters->pi/lambda_perturb*(*grid)(i,j,k).y);
          (*(*phi)(1))(i,j,k) = 1-.5*ratio*tanh(2.*((*grid)(i,j,k).z - zeta + 
                parameters->z_length/interface_loc)/delta*atanh(alpha));          
        }
      }
    }
  }
  if(parameters->progressive_wave){
    //Progressive wave 
    //density
    for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++) {
      for(int j=grid->J_Min_With_Halo(); j<=grid->J_Max_With_Halo(); j++) {
        for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++) {
          zeta = delta_perturb*cos(2*parameters->pi/lambda_perturb*(*grid)(i,j,k).y);
          (*(*phi)(1))(i,j,k) = 1-.5*ratio*tanh(2.*((*grid)(i,j,k).z - zeta + 
                parameters->z_length/interface_loc)/delta*atanh(alpha));           
        }
      }
    }
  }
  //passive scalar
  if(parameters->num_scalars == 2){
    for(int i=grid->I_Min_With_Halo(); i<=grid->I_Max_With_Halo(); i++) {
      for(int j=grid->J_Min_With_Halo(); j<=grid->J_Max_With_Halo(); j++) {
        for(int k=grid->K_Min_With_Halo(); k<=grid->K_Max_With_Halo(); k++) {
          if(((*grid)(i,j,k).x > 2.25) && ((*grid)(i,j,k).x < 2.75))
            (*(*phi)(2))(i,j,k) = 1.;
        }
      }
    }
  }
}

//POSTPROCESS: enforce BCs and populate ghost cells
assert(mpi_driver);
Set_Progressive_Wave_BC(parameters->time);
Enforce_Velocity_BC(*u);  
mpi_driver->Exchange_Ghost_Values_For_Vector_Field(*u); 

/*
//test print
if(mpi_driver->west_proc == MPI_PROC_NULL) {
    cout << "Time step 0" << endl;
    for(int j = grid->J_Max(); j >= grid->J_Min(); j--) {
      //cout << ( (*u)(grid->I_Min(),j,ceil(grid->K_Max()/2)).x + 
      //          (*u)(grid->I_Min()-1,j,ceil(grid->K_Max()/2)).x ) / 2. << endl;
      cout << (*u)(grid->I_Min()-2,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()-1,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()  ,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()+1,j,ceil(grid->K_Max()/2)).x << "\t" <<
              (*u)(grid->I_Min()+2,j,ceil(grid->K_Max()/2)).x << "\t" <<
              ( (*u)(grid->I_Min(),j,ceil(grid->K_Max()/2)).x + 
                (*u)(grid->I_Min()-1,j,ceil(grid->K_Max()/2)).x ) / 2. << endl;
    }
  }
*/
}
//*****************************************************************************
// Set progressive wave BC
//*****************************************************************************
  template<class T> 
void NAVIER_STOKES_SOLVER<T>::Set_Progressive_Wave_BC(T time)
{
  if(parameters->west_velocity)
    for(int j=parameters->j_min_w_h; j<=parameters->j_max_w_h; j++)
      for(int k=parameters->k_min_w_h; k<=parameters->k_max_w_h; k++) 
        (*parameters->west_velocity)(j,k).x = 
             parameters->forcing_amp*cos(parameters->m*(*grid)(grid->I_Min(),j,k).z)
                                    *sin(parameters->freq*time);
} 
//*****************************************************************************
// Total Kinetic Energy
//*****************************************************************************
  template<class T> 
T NAVIER_STOKES_SOLVER<T>::Total_Kinetic_Energy()
{
  T E_k = (T)0;
  T F_Ek = (T)0;
  // calculate local E_k and F_Ek
  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
        T cell_volume = (T)1 / (*grid->inverse_Jacobian)(i,j,k);
        E_k += 0.5 * cell_volume * 
               ( pow((*u)(i,j,k).x,2) + pow((*u)(i,j,k).y,2) + pow((*u)(i,j,k).z,2) ); 
      }
  if(parameters->west_velocity)
    if(mpi_driver->west_proc == NULL)
      for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
        for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
          int i = grid->I_Min();
          T cell_volume = (T)1 / (*grid->inverse_Jacobian)(i,j,k);
          T cell_area = cell_volume / (parameters->x_length/parameters->num_total_nodes_x);
          F_Ek += 0.5 * cell_area * (*parameters->west_velocity)(j,k).x *
            ( pow((*parameters->west_velocity)(j,k).x,2) 
              + pow((*u)(i,j,k).y,2) + pow((*u)(i,j,k).z,2) ); 
         }
  // sum over all procs
  mpi_driver->Replace_With_Sum_On_All_Procs(E_k);
  mpi_driver->Replace_With_Sum_On_All_Procs(F_Ek);
  
  //write to disk
  if(!mpi_driver->my_rank){

    stringstream filename;

    //overwrite existing file
    if(parameters->time_step==parameters->save_data_timestep_period){
      filename << parameters->output_dir << "kinetic_energy.0";
      ofstream output(filename.str().c_str(), ios::out | ios::trunc | ios::binary);
      if(!output){
        cout << "ERROR: could not open potential_energy file for writing" <<endl;
        return 0;
      }
      output.write(reinterpret_cast<char *>(&E_k),sizeof(T));
      output.write(reinterpret_cast<char *>(&F_Ek),sizeof(T));
      output.close();
    } 

    //append at end of existing file
    else{
      filename << parameters->output_dir << "kinetic_energy.0";
      ofstream output(filename.str().c_str(), ios::out | ios::app | ios::binary);
      if(!output){
        cout << "ERROR: could not open potential_energy file for writing" <<endl;
        return 0;
      }
      output.write(reinterpret_cast<char *>(&E_k),sizeof(T));
      output.write(reinterpret_cast<char *>(&F_Ek),sizeof(T));
      output.close();
    }
  }
  return E_k;
}
//*****************************************************************************
template class NAVIER_STOKES_SOLVER<double>;
