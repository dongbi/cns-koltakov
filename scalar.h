// Scalar field evolution class

#ifndef __SCALAR__
#define __SCALAR__

template<class T=double>
class SCALAR
{
 public:
  SCALAR(PARAMETERS<T> *par, MPI_DRIVER<T> *md, CURVILINEAR_GRID<T> *g, 
	 CONVECTION<T> *c, TURBULENCE<T> *t, ARRAY_3D<T> *rho);
  ~SCALAR();

  void Update_RHS();
  void Solve();

  T Rho_Rest(const int i, const int j, const int k){return (*Rho_rest)(i,j,k);}
  void Set_Initial_Density_Profile();
  void Set_Uniform_Density_Profile(const T rho_const) 
   {Rho_rest->Set_All_Elements_To(rho_const); *Rho = *Rho_rest;}

 private:
  void Enforce_Density_BC(ARRAY_3D<T>& rho);

  PARAMETERS<T>* parameters;
  MPI_DRIVER<T>* mpi_driver;
  CURVILINEAR_GRID<T>* grid;
  CONVECTION<T>* convection;
  TURBULENCE<T>* turbulence;
  
  ARRAY_3D<T> *Rho, *Rho_rest, *RHS, *RHS_for_AB;
  T imin, imax, jmin, jmax, kmin, kmax, halo; //set from grid
};
//*****************************************************************************
// Constructor
//*****************************************************************************
template<class T>
SCALAR<T>::SCALAR(PARAMETERS<T> *par,MPI_DRIVER<T> *md,CURVILINEAR_GRID<T> *g,
		  CONVECTION<T> *c, TURBULENCE<T> *t, ARRAY_3D<T> *rho)
: parameters(par), mpi_driver(md), grid(g), convection(c),turbulence(t),Rho(rho)
{
  imin = g->I_Min(); imax = g->I_Max(); jmin = g->J_Min(); jmax = g->J_Max();
  kmin = g->K_Min(); kmax = g->K_Max(); halo = g->Halo_Size();

  RHS        = new ARRAY_3D<T>(imin,imax, jmin,jmax, kmin,kmax, halo);
  RHS_for_AB = new ARRAY_3D<T>(*RHS);
  Rho_rest = new ARRAY_3D<T>(*rho);

  //Set_Initial_Density_Profile(); // setup a sloshing wave

  //Set_Uniform_Density_Profile((T)1);
}
//*****************************************************************************
// Destructor
//*****************************************************************************
template<class T>
SCALAR<T>::~SCALAR()
{
  delete RHS; delete RHS_for_AB; delete Rho_rest;
}
//*****************************************************************************
// Resets the RHS for the scalar (rho) evolution linear system
//*****************************************************************************
template<class T> 
void SCALAR<T>::Update_RHS()
{
  T leading_term_AB_coefficient;

  // Forward Euler timesteping on the first time step
  if(parameters->time_step == 1)
    leading_term_AB_coefficient = (T)1;
  else{
    leading_term_AB_coefficient = (T)1.5;
    // adding Adams-Bashforth contribution from previous time step
    *RHS = (T)-.5 * (*RHS_for_AB);
  }
  // reset array for the new time step
  RHS_for_AB->Set_All_Elements_To((T)0);

  // convection term with SHARP
  convection->Add_Scalar_Convection_Term(*RHS_for_AB);
  //mpi_driver->Write_Global_Array_To_Disk("convection_scalar_rhs",*RHS_for_AB, 
  //					 parameters->time_step);
  // add convection term to correct for moving grid velocity
  if(parameters->moving_grid)
    convection->Add_Moving_Grid_Convection_Term(*RHS_for_AB);  

  // LES self-similarity term
  if(parameters->turbulence){
    //RHS_for_AB += *turbulence->tau; //need 4th comp
  }

  // Cross diffusive terms at time = n-1 discretized with Crank-Nicolson
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++){
	//D_12 and D_13
	T diff_minus   = parameters->molecular_diffusivity,
	  diff_plus    = parameters->molecular_diffusivity,
	  diff_term;
        if(parameters->turbulence){
	  diff_minus += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
			        (*turbulence->eddy_diffusivity)(i-1,j,k) );
	  diff_plus  += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
	                        (*turbulence->eddy_diffusivity)(i+1,j,k) );
	}
	diff_term  = diff_plus * (*grid->G12)(i  ,j,k)*
                                      ( (*Rho)(i  ,j+1,k) - (*Rho)(i  ,j-1,k) 
                                      + (*Rho)(i+1,j+1,k) - (*Rho)(i+1,j-1,k) );
	          - diff_minus* (*grid->G12)(i-1,j,k)*
	                              ( (*Rho)(i  ,j+1,k) - (*Rho)(i  ,j-1,k) 
                                      + (*Rho)(i-1,j+1,k) - (*Rho)(i-1,j-1,k) );
	diff_term += diff_plus * (*grid->G13)(i  ,j,k)*
	                              ( (*Rho)(i  ,j,k+1) - (*Rho)(i  ,j,k-1) 
                                      + (*Rho)(i+1,j,k+1) - (*Rho)(i+1,j,k-1) );
	          - diff_minus* (*grid->G13)(i-1,j,k)*
	                              ( (*Rho)(i  ,j,k+1) - (*Rho)(i  ,j,k+1) 
                                      + (*Rho)(i-1,j,k+1) - (*Rho)(i-1,j,k-1) );

	//D_23 and D_21
	diff_minus = parameters->molecular_diffusivity;
	diff_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diff_minus += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
			        (*turbulence->eddy_diffusivity)(i,j-1,k) );
	  diff_plus  += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
	                        (*turbulence->eddy_diffusivity)(i,j+1,k) );
	}
	diff_term += diff_plus * (*grid->G23)(i,j  ,k)*
	                              ( (*Rho)(i,j  ,k+1) - (*Rho)(i,j  ,k-1) 
                                      + (*Rho)(i,j+1,k+1) - (*Rho)(i,j+1,k-1) );
	           - diff_minus* (*grid->G23)(i,j-1,k)*
	                              ( (*Rho)(i,j  ,k+1) - (*Rho)(i,j  ,k-1) 
                                      + (*Rho)(i,j-1,k+1) - (*Rho)(i,j-1,k-1) );
	diff_term += diff_plus * (*grid->G21)(i,j  ,k)*
	                              ( (*Rho)(i+1,j  ,k) - (*Rho)(i-1,j  ,k) 
                                      + (*Rho)(i+1,j+1,k) - (*Rho)(i-1,j+1,k) );
	           - diff_minus* (*grid->G21)(i,j-1,k)*
	                              ( (*Rho)(i+1,j  ,k) - (*Rho)(i-1,j  ,k) 
                                      + (*Rho)(i+1,j-1,k) - (*Rho)(i-1,j-1,k) );

	//D_31 and D_32
	diff_minus = parameters->molecular_diffusivity;
	diff_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diff_minus += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
			        (*turbulence->eddy_diffusivity)(i,j,k-1) );
	  diff_plus  += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
	                        (*turbulence->eddy_diffusivity)(i,j,k+1) );
	}
	diff_term += diff_plus * (*grid->G31)(i,j,k  )*
	                              ( (*Rho)(i+1,j,k  ) - (*Rho)(i-1,j,k  ) 
                                      + (*Rho)(i+1,j,k+1) - (*Rho)(i-1,j,k+1) );
	           - diff_minus* (*grid->G31)(i,j,k-1)*
	                              ( (*Rho)(i+1,j,k  ) - (*Rho)(i-1,j,k  ) 
                                      + (*Rho)(i+1,j,k-1) - (*Rho)(i-1,j,k-1) );
	diff_term += diff_plus * (*grid->G32)(i,j,k  )*
	                              ( (*Rho)(i,j+1,k  ) - (*Rho)(i,j-1,k  ) 
                                      + (*Rho)(i,j+1,k+1) - (*Rho)(i,j-1,k+1) );
	           - diff_minus* (*grid->G32)(i,j,k-1)*
	                              ( (*Rho)(i,j+1,k  ) - (*Rho)(i,j-1,k  ) 
                                      + (*Rho)(i,j+1,k-1) - (*Rho)(i,j-1,k-1) );

	// add D_12+D_13+D_23+D_21+D_31+D_32 to RHS
        (*RHS_for_AB)(i,j,k) += diff_term;
  }
  
  // adding Adams-Bashforth contribution for current time step
  // 'RHS_for_AB' is saved from this point until the next time step
  *RHS += leading_term_AB_coefficient * (*RHS_for_AB);

  // Diagonal viscous terms at time = n-1 discretized with Crank-Nicolson
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++){
	//D_11
	T diff_minus = parameters->molecular_diffusivity,
	  diff_plus  = parameters->molecular_diffusivity,
	  diff_term;
        if(parameters->turbulence){
	  diff_minus += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
			        (*turbulence->eddy_diffusivity)(i-1,j,k) );
	  diff_plus  += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
	                        (*turbulence->eddy_diffusivity)(i+1,j,k) );
	}
	diff_term = 
          diff_plus * (*grid->G11)(i  ,j,k) * ((*Rho)(i+1,j,k)-(*Rho)(i  ,j,k))
	- diff_minus* (*grid->G11)(i-1,j,k) * ((*Rho)(i  ,j,k)-(*Rho)(i-1,j,k));
	//D_22
	diff_minus = parameters->molecular_diffusivity;
	diff_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diff_minus += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
			        (*turbulence->eddy_diffusivity)(i,j-1,k) );
	  diff_plus  += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
	                        (*turbulence->eddy_diffusivity)(i,j+1,k) );
	}
	diff_term +=
          diff_plus * (*grid->G22)(i,j  ,k) * ((*Rho)(i,j+1,k)-(*Rho)(i,j  ,k))
	- diff_minus* (*grid->G22)(i,j-1,k) * ((*Rho)(i,j  ,k)-(*Rho)(i,j-1,k));
	//D_33
	diff_minus = parameters->molecular_diffusivity;
	diff_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diff_minus += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
			        (*turbulence->eddy_diffusivity)(i,j,k-1) );
	  diff_plus  += (T).5*( (*turbulence->eddy_diffusivity)(i,j,k) + 
	                        (*turbulence->eddy_diffusivity)(i,j,k+1) );
	}
	diff_term +=
          diff_plus * (*grid->G33)(i,j,k  ) * ((*Rho)(i,j,k+1)-(*Rho)(i,j,k  ))
	- diff_minus* (*grid->G33)(i,j,k-1) * ((*Rho)(i,j,k  )-(*Rho)(i,j,k-1));
	// add D_11+D_22+D_33 to RHS
	(*RHS)(i,j,k) += diff_term;
  }

  // Jacobian change in one time_step adjustment term for a moving grid
  if(parameters->moving_grid)
    *RHS -= (*((CURVILINEAR_MOVING_GRID<T>*)grid)->Jacobian_diff) * (*Rho);

  // Multiply every term by delta_time and inverse Jacobian
  *RHS *=  parameters->delta_time * (*grid->inverse_Jacobian); 
  //mpi_driver->Write_Global_Array_To_Disk("presolve_scalar_rhs", *RHS, 
  //					 parameters->time_step);
}
//*****************************************************************************
// Solves scalar (rho) evolution linear system using member var RHS for rhs
//*****************************************************************************
template<class T>
void SCALAR<T>::Solve()
{
  TRIDIAGONAL_SOLVER<T> LS_Solver(*mpi_driver);
  ARRAY_2D<T> *A, *B, *C, *F;  // LHS diagonals & RHS for tridiagonal solve

  // Solve 3 linear systems for scalar evolution equation: M_i*M_j*M_k*x = RHS
  //--------------------------------------------------------------------------
  // Solve for I-DIRECTION: (I-(dt/2J^{-1})*D_11)(rho^*-rho^n)_i = RHS
  //--------------------------------------------------------------------------
  A = new ARRAY_2D<T>(jmin,jmax, imin-1,imax+1, 0); //halo=0
  B = new ARRAY_2D<T>(jmin,jmax, imin-1,imax+1, 0);
  C = new ARRAY_2D<T>(jmin,jmax, imin-1,imax+1, 0);
  F = new ARRAY_2D<T>(jmin,jmax, imin-1,imax+1, 0);

  for(int k = kmin; k <= kmax; k++){
    // Construct tridiagonal matrices A,B,C for LHS and F for RHS(of system)
    for(int j = jmin; j <= jmax; j++)
      for(int i = imin; i <= imax; i++){
	T diffusivity_minus = parameters->molecular_diffusivity,
	  diffusivity_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diffusivity_minus += (T).5*((*turbulence->eddy_diffusivity)(i,j,k) + 
			              (*turbulence->eddy_diffusivity)(i-1,j,k));
	  diffusivity_plus  += (T).5*((*turbulence->eddy_diffusivity)(i,j,k) + 
	                              (*turbulence->eddy_diffusivity)(i+1,j,k));
	}
	(*A)(j,i) = - diffusivity_minus * (T).5 * parameters->delta_time * 
	            (*grid->inverse_Jacobian)(i,j,k) * (*grid->G11)(i-1,j,k);
	(*C)(j,i) = - diffusivity_plus  * (T).5 * parameters->delta_time * 
	            (*grid->inverse_Jacobian)(i,j,k) * (*grid->G11)(i,j,k);
        (*B)(j,i) = (T)1 - (*A)(j,i) - (*C)(j,i);
	(*F)(j,i) = (*RHS)(i,j,k);
      }// for: j,i
    // Boundary conditions for I-direction: west
    if(mpi_driver->west_proc == MPI_PROC_NULL)
      for(int j = jmin; j <= jmax; j++){
	(*A)(j,imin-1) = (T)0;
	(*B)(j,imin-1) = (T)1;
	(*C)(j,imin-1) = (T)-1;
	// RHS : BC west
	(*RHS_for_AB)(imin-1,j,k) = (*grid->G12)(imin-1,j,k) * 
                                ( (*Rho)(imin-1,j+1,k) - (*Rho)(imin-1,j-1,k)
                                + (*Rho)(imin  ,j+1,k) - (*Rho)(imin  ,j-1,k) )
	                          + (*grid->G13)(imin-1,j,k) * 
                                ( (*Rho)(imin-1,j,k+1) - (*Rho)(imin-1,j,k-1) 
                                + (*Rho)(imin  ,j,k+1) - (*Rho)(imin  ,j,k-1) );
	(*RHS_for_AB)(imin-1,j,k) /= (*grid->G11)(imin-1,j,k);
	(*F)(j,imin-1) = (*RHS_for_AB)(imin-1,j,k);
    }// west BC
    // Boundary conditions for I-direction: east
    if(mpi_driver->east_proc == MPI_PROC_NULL)
      for(int j = jmin; j <= jmax; j++){
	(*A)(j,imax+1) = (T)1;
	(*B)(j,imax+1) = (T)-1;
	(*C)(j,imax+1) = (T)0;
	// RHS : BC east
	(*RHS_for_AB)(imax+1,j,k) = (*grid->G12)(imax,j,k) * 
                                ( (*Rho)(imax  ,j+1,k) - (*Rho)(imax  ,j-1,k)
                                + (*Rho)(imax+1,j+1,k) - (*Rho)(imax+1,j-1,k) )
	                          + (*grid->G13)(imax,j,k) * 
                                ( (*Rho)(imax  ,j,k+1) - (*Rho)(imax  ,j,k-1) 
                                + (*Rho)(imax+1,j,k+1) - (*Rho)(imax+1,j,k-1) );
	(*RHS_for_AB)(imax+1,j,k) /= (*grid->G11)(imax,j,k);
	(*F)(j,imax+1) = (*RHS_for_AB)(imax+1,j,k);
    }//east BC

    //solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(*A, *B, *C, *F, 
      parameters->periodic_in_x, mpi_driver->west_proc, mpi_driver->east_proc);

    //solution of linear system is in F
    for(int j = jmin; j <= jmax; j++)
      for(int i = imin; i <= imax; i++)
	(*RHS)(i,j,k) = (*F)(j,i);
  }//for: k
  delete A; delete B; delete C; delete F;

  //--------------------------------------------------------------------------
  // Solve for J-DIRECTION: (I-(dt/2J^{-1})*D_22)(rho^*-rho^n)_j = RHS 
  //--------------------------------------------------------------------------
  A = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0); //halo=0
  B = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0);
  C = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0);
  F = new ARRAY_2D<T>(imin,imax, jmin-1,jmax+1, 0);

  for(int k = kmin; k <= kmax; k++){
    // Construct tridiagonal matrices A,B,C for LHS and F for RHS (of LinSys)
    for(int i = imin; i <= imax; i++)
      for(int j = jmin; j <= jmax; j++){
	T diffusivity_minus = parameters->molecular_diffusivity,
	  diffusivity_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diffusivity_minus += (T).5*((*turbulence->eddy_diffusivity)(i,j  ,k)
                                   +  (*turbulence->eddy_diffusivity)(i,j-1,k));
	  diffusivity_plus  += (T).5*((*turbulence->eddy_diffusivity)(i,j  ,k) 
                                    + (*turbulence->eddy_diffusivity)(i,j+1,k));
	}
	(*A)(i,j) = - diffusivity_minus * (T).5 * parameters->delta_time * 
	            (*grid->inverse_Jacobian)(i,j,k) * (*grid->G22)(i,j-1,k);
	(*C)(i,j) = - diffusivity_plus  * (T).5 * parameters->delta_time * 
	            (*grid->inverse_Jacobian)(i,j,k) * (*grid->G22)(i,j  ,k);
        (*B)(i,j) = (T)1 - (*A)(i,j) - (*C)(i,j);
	(*F)(i,j) = (*RHS)(i,j,k);
      }// for: i,j
    // Boundary conditions for J-direction: south
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
	(*A)(i,jmin-1) = (T)0;
	(*B)(i,jmin-1) = (T)1;
	(*C)(i,jmin-1) = (T)-1;//0
	// RHS : south BC
	(*RHS_for_AB)(i,jmin-1,k) = (*grid->G23)(i,jmin-1,k) * 
	                        ( (*Rho)(i,jmin-1,k+1) - (*Rho)(i,jmin-1,k-1)
				+ (*Rho)(i,jmin  ,k+1) - (*Rho)(i,jmin  ,k-1) )
	                          + (*grid->G21)(i,jmin-1,k) * 
	                        ( (*Rho)(i+1,jmin-1,k) - (*Rho)(i-1,jmin-1,k) 
				+ (*Rho)(i+1,jmin  ,k) - (*Rho)(i-1,jmin  ,k) );
	(*RHS_for_AB)(i,jmin-1,k) /= (*grid->G22)(i,jmin-1,k);
	(*F)(i,jmin-1) = (*RHS_for_AB)(i,jmin-1,k);
    }
    // Boundary conditions for J-direction: north
    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
	(*A)(i,jmax+1) = (T)1;//0;	
	(*B)(i,jmax+1) = (T)-1;//1;
	(*C)(i,jmax+1) = (T)0;
	// RHS : north BC
	(*RHS_for_AB)(i,jmax+1,k) = (*grid->G23)(i,jmax,k) * 
	                        ( (*Rho)(i,jmax  ,k+1) - (*Rho)(i,jmax  ,k-1)
				+ (*Rho)(i,jmax+1,k+1) - (*Rho)(i,jmax+1,k-1) )
	                          + (*grid->G21)(i,jmax,k) * 
	                        ( (*Rho)(i+1,jmax  ,k) - (*Rho)(i-1,jmax  ,k) 
				+ (*Rho)(i+1,jmax+1,k) - (*Rho)(i-1,jmax+1,k) );
	(*RHS_for_AB)(i,jmax+1,k) /= (*grid->G22)(i,jmax,k);
	(*F)(i,jmax+1) = (*RHS_for_AB)(i,jmax+1,k);
    }
 
    //solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(*A, *B, *C, *F, 
      parameters->periodic_in_y, mpi_driver->suth_proc, mpi_driver->nrth_proc);

    //solution of linear system is in F
    for(int i = imin; i <= imax; i++)
      for(int j = jmin; j <= jmax; j++)
	(*RHS)(i,j,k) = (*F)(i,j);
  }//for: k
  delete A; delete B; delete C; delete F;  

  //--------------------------------------------------------------------------
  // Solve for K-DIRECTION: (I-(dt/2J^{-1})*D_33)(rho^*-rho^n)_k = RHS
  //--------------------------------------------------------------------------
  A = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0); //halo=0
  B = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0);
  C = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0);
  F = new ARRAY_2D<T>(imin,imax, kmin-1,kmax+1, 0);

  for(int j = jmin; j <= jmax; j++){
    // Construct tridiagonal matrices A,B,C for LHS and F for RHS (of LinSys)
    for(int k = kmin; k <= kmax; k++)
      for(int i = imin; i <= imax; i++){
	T diffusivity_minus = parameters->molecular_diffusivity,
	  diffusivity_plus  = parameters->molecular_diffusivity;
        if(parameters->turbulence){
	  diffusivity_minus += (T).5*((*turbulence->eddy_diffusivity)(i,j,k  )  
			             +(*turbulence->eddy_diffusivity)(i,j,k-1));
	  diffusivity_plus  += (T).5*((*turbulence->eddy_diffusivity)(i,j,k  ) 
	                             +(*turbulence->eddy_diffusivity)(i,j,k+1));
	}
	(*A)(i,k) = - diffusivity_minus * (T).5 * parameters->delta_time * 
	            (*grid->inverse_Jacobian)(i,j,k) * (*grid->G33)(i,j,k-1);
	(*C)(i,k) = - diffusivity_plus  * (T).5 * parameters->delta_time * 
	            (*grid->inverse_Jacobian)(i,j,k) * (*grid->G33)(i,j,k  );
        (*B)(i,k) = (T)1 - (*A)(i,k) - (*C)(i,k);
	(*F)(i,k) = (*RHS)(i,j,k);
    }// for: i,k
    // Boundary conditions for K-direction: back
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
	(*A)(i,kmin-1) = (T)0;
	(*B)(i,kmin-1) = (T)1;
	(*C)(i,kmin-1) = (T)-1; 
	// RHS : BC back
	(*RHS_for_AB)(i,j,kmin-1) = (*grid->G31)(i,j,kmin-1) * 
	                        ( (*Rho)(i+1,j,kmin-1) - (*Rho)(i-1,j,kmin-1)
	                        + (*Rho)(i+1,j,kmin  ) - (*Rho)(i-1,j,kmin  ) )
	                          + (*grid->G32)(i,j,kmin-1) * 
	                        ( (*Rho)(i,j+1,kmin-1) - (*Rho)(i,j-1,kmin-1) 
	                        + (*Rho)(i,j+1,kmin  ) - (*Rho)(i,j-1,kmin  ) );
	(*RHS_for_AB)(i,j,kmin-1) /= (*grid->G33)(i,j,kmin-1);
	(*F)(i,kmin-1) = (*RHS_for_AB)(i,j,kmin-1);
    }// back BC
    // Boundary conditions for K-direction: front
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int i = imin; i <= imax; i++){
	(*A)(i,kmax+1) = (T)1;
	(*B)(i,kmax+1) = (T)-1;
	(*C)(i,kmax+1) = (T)0;
	// RHS : BC front
	(*RHS_for_AB)(i,j,kmax+1) = (*grid->G31)(i,j,kmax) * 
	                        ( (*Rho)(i+1,j,kmax  ) - (*Rho)(i-1,j,kmax  )
	                        + (*Rho)(i+1,j,kmax+1) - (*Rho)(i-1,j,kmax+1) )
	                          + (*grid->G32)(i,j,kmax) * 
	                        ( (*Rho)(i,j+1,kmax  ) - (*Rho)(i,j-1,kmax  ) 
	                        + (*Rho)(i,j+1,kmax+1) - (*Rho)(i,j-1,kmax+1) );
	(*RHS_for_AB)(i,j,kmax+1) /= (*grid->G33)(i,j,kmax);
	(*F)(i,kmax+1) = (*RHS_for_AB)(i,j,kmax+1);
    }//front BC

    //solve tridiagonal system
    LS_Solver.Solve_Array_Of_Tridiagonal_Linear_Systems(*A, *B, *C, *F, 
      parameters->periodic_in_z, mpi_driver->back_proc, mpi_driver->frnt_proc);

    //solution of linear system is in F
    for(int i = imin; i <= imax; i++)
      for(int k = kmin; k <= kmax; k++)
	(*RHS)(i,j,k) = (*F)(i,k);
  }//for: k
  delete A; delete B; delete C; delete F;

  // Update scalar field
  //*Rho += *RHS;
  for(int i = imin; i <= imax; i++)
    for(int j = jmin; j <= jmax; j++)
      for(int k = kmin; k <= kmax; k++)
	(*Rho)(i,j,k) += (*RHS)(i,j,k);

  //mpi_driver->Write_Global_Array_To_Disk("scalar_rhs", *RHS, 
  // 					 parameters->time_step);

  //ut<<"RHS_center="<<(*RHS)(imax,round(.5*jmax),round(.5*kmax))<<endl;
  //ut<<"rho_halo="<<(*Rho)(imax+2,round(.5*jmax),round(.5*kmax))<<endl;

  // BCs for Rho
  Enforce_Density_BC(*Rho);
 
  // Swap halo regions among procs
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*Rho);  
 
  mpi_driver->Syncronize_All_Procs(); // Wait for all procs to finish
}
//*****************************************************************************
template<class T>
void SCALAR<T>::Enforce_Density_BC(ARRAY_3D<T>& rho)
{
  // face BC: Left
  if(mpi_driver->west_proc == MPI_PROC_NULL)
    for(int j=jmin-halo; j<=jmax+halo; j++)
      for(int k=kmin-halo; k<=kmax+halo; k++){
	rho(imin-1,j,k) = rho(imin,j,k) + (*RHS_for_AB)(imin-1,j,k);
	rho(imin-2,j,k) = (T)3 * (rho(imin-1,j,k) - rho(imin,j,k))
                        +  rho(imin+1,j,k);
  }
  // face BC: Right
  if(mpi_driver->east_proc == MPI_PROC_NULL)
    for(int j=jmin-halo; j<=jmax+halo; j++)
      for(int k=kmin-halo; k<=kmax+halo; k++){
	rho(imax+1,j,k) = rho(imax,j,k) - (*RHS_for_AB)(imax+1,j,k);
	rho(imax+2,j,k) = (T)3 * (rho(imax+1,j,k) - rho(imax,j,k))
                        +  rho(imax-1,j,k);
  }
  // face BC: Bottom
  if(mpi_driver->suth_proc == MPI_PROC_NULL)
    for(int i=imin-halo; i<=imax+halo; i++)
      for(int k=kmin-halo; k<=kmax+halo; k++){
	rho(i,jmin-1,k) = rho(i,jmin,k) + (*RHS_for_AB)(i,jmin-1,k);
	rho(i,jmin-2,k) = (T)3 * (rho(i,jmin-1,k) - rho(i,jmin,k))
	                +  rho(i,jmin+1,k);
	//rho(i,jmin-2,k) = rho(i,jmin-1,k) = rho(i,jmin,k); /////
  }
  // face BC: Top
  if(mpi_driver->nrth_proc == MPI_PROC_NULL)
    for(int i=imin-halo; i<=imax+halo; i++)
      for(int k=kmin-halo; k<=kmax+halo; k++){
	rho(i,jmax+1,k) = rho(i,jmax,k) - (*RHS_for_AB)(i,jmax+1,k);
	rho(i,jmax+2,k) = (T)3 * (rho(i,jmax+1,k) - rho(i,jmax,k))
	                +  rho(i,jmax-1,k);
  }
  // face BC: Back
  if(mpi_driver->back_proc == MPI_PROC_NULL)
    for(int i=imin-halo; i<=imax+halo; i++)
      for(int j=jmin-halo; j<=jmax+halo; j++){
	rho(i,j,kmin-1) = rho(i,j,kmin) + (*RHS_for_AB)(i,j,kmin-1);
	rho(i,j,kmin-2) = (T)3 * (rho(i,j,kmin-1) - rho(i,j,kmin))
	                +  rho(i,j,kmin+1);
  }
  // face BC: Front
  if(mpi_driver->frnt_proc == MPI_PROC_NULL)
    for(int i=imin-halo; i<=imax+halo; i++)
      for(int j=jmin-halo; j<=jmax+halo; j++){
	rho(i,j,kmax+1) = rho(i,j,kmax) - (*RHS_for_AB)(i,j,kmax+1);
	rho(i,j,kmax+2) = (T)3 * (rho(i,j,kmax+1) - rho(i,j,kmax))
	                +  rho(i,j,kmax-1);
  }
  // edge BC: 4 edges on left face
  if(mpi_driver->west_proc == MPI_PROC_NULL) {
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int k=kmin; k<=kmax; k++)
	rho(imin-1,jmin-1,k) = (T).5 * (rho(imin,jmin-1,k)+rho(imin-1,jmin,k));
    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int k=kmin; k<=kmax; k++)
	rho(imin-1,jmax+1,k) = (T).5 * (rho(imin,jmax+1,k)+rho(imin-1,jmax,k));
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int j=jmin; j<=jmax; j++)
	rho(imin-1,j,kmin-1) = (T).5 * (rho(imin,j,kmin-1)+rho(imin-1,j,kmin));
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int j=jmin; j<=jmax; j++)
	rho(imin-1,j,kmax+1) = (T).5 * (rho(imin,j,kmax+1)+rho(imin-1,j,kmax));
  }
  // edge BC: 4 edges on right face
  if(mpi_driver->east_proc == MPI_PROC_NULL) {
    if(mpi_driver->suth_proc == MPI_PROC_NULL)
      for(int k=kmin; k<=kmax; k++)
	rho(imax+1,jmin-1,k) = (T).5 * (rho(imax,jmin-1,k)+rho(imax+1,jmin,k));
    if(mpi_driver->nrth_proc == MPI_PROC_NULL)
      for(int k=kmin; k<=kmax; k++)
	rho(imax+1,jmax+1,k) = (T).5 * (rho(imax,jmax+1,k)+rho(imax+1,jmax,k));
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int j=jmin; j<=jmax; j++)
	rho(imax+1,j,kmin-1) = (T).5 * (rho(imax,j,kmin-1)+rho(imax+1,j,kmin));
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int j=jmin; j<=jmax; j++)
	rho(imax+1,j,kmax+1) = (T).5 * (rho(imax,j,kmax+1)+rho(imax+1,j,kmax));
  }
  // edge BC: 2 streamwise edges on bottom face
  if(mpi_driver->suth_proc == MPI_PROC_NULL) {
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int i=imin; i<=imax; i++)
	rho(i,jmin-1,kmin-1) = (T).5 * (rho(i,jmin,kmin-1)+rho(i,jmin-1,kmin));
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int i=imin; i<=imax; i++)
	rho(i,jmin-1,kmax+1) = (T).5 * (rho(i,jmin,kmax+1)+rho(i,jmin-1,kmax));
  }
  // edge BC: 2 streamwise edges on top face
  if(mpi_driver->nrth_proc == MPI_PROC_NULL) {
    if(mpi_driver->back_proc == MPI_PROC_NULL)
      for(int i=imin; i<=imax; i++)
	rho(i,jmax+1,kmin-1) = (T).5 * (rho(i,jmax,kmin-1)+rho(i,jmax+1,kmin));
    if(mpi_driver->frnt_proc == MPI_PROC_NULL)
      for(int i=imin; i<=imax; i++)
	rho(i,jmax+1,kmax+1) = (T).5 * (rho(i,jmax,kmax+1)+rho(i,jmax+1,kmax));
  }
}
//*****************************************************************************
// Initializes the density (Rho) array
//*****************************************************************************
template<class T>
void SCALAR<T>::Set_Initial_Density_Profile()
{
  T kk = parameters->pi, ka = .1, kd = parameters->pi, alpha = .99, 
    k_delta = .05*parameters->pi, rho0 = 1000., delta_rho = .03*rho0;

  for(int i = imin-halo; i <= imax+halo; i++)
    for(int j = jmin-halo; j <= jmax+halo; j++)
      for(int k = kmin-halo; k <= kmax+halo; k++){
	T k_eta = ka * ( (1.-ka*ka/64.)*cos(kk*(*grid)(i,j,k).x) 
			    -(ka*ka/8.)*cos(3.*kk*(*grid)(i,j,k).x) );
	(*Rho)(i,j,k) =-.5*delta_rho* tanh(2.* atanh(alpha)
                                             * (kk*(*grid)(i,j,k).y-k_eta+kd/2.)
                                             / k_delta);
      }
  /*
  int half_j = round((jmax-jmin+1)/2);
  for(int i = imin-halo; i <= imax+halo; i++)
    for(int j = jmin-halo; j <= jmax+halo; j++)
      for(int k = kmin-halo; k <= kmax+halo; k++)
	if(j<half_j) (*Rho_rest)(i,j,k) = (*Rho)(imin,jmin,kmin); 
	else (*Rho_rest)(i,j,k) = (*Rho)(imax,jmax,kmax); 
  */
}
//*****************************************************************************
#endif
