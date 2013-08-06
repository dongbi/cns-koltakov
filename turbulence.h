// turbulence.h

#ifndef __TURBULENCE__
#define __TURBULENCE__

#include "array_3d.h"

template <class T=double> 
class TURBULENCE
{
 public:
  TURBULENCE(PARAMETERS<T> *par, MPI_DRIVER<T> *md, CURVILINEAR_GRID<T> *g,
      ARRAY_3D<T> *rhoin, ARRAY_3D<VECTOR_3D<T> > *uin);
  ~TURBULENCE();

  void Calculate_Eddy_Viscosity();

  ARRAY_3D<T> *eddy_viscosity, *eddy_diffusivity;
  ARRAY_3D<VECTOR_3D<T> > *tau;

 private:
  PARAMETERS<T> *parameters;
  MPI_DRIVER<T> *mpi_driver;
  CURVILINEAR_GRID<T> *grid;
  ARRAY_3D<T> *rho;
  ARRAY_3D<VECTOR_3D<T> > *u;

  int i, j, k, m;
  int imin, imax, jmin, jmax, kmin, kmax, halo;
  T tiny, alpha2, weit, fact, temp, fmin, fmax;
  ARRAY_1D<ARRAY_3D<T> > *ss, *tt, *zz, *rr;

};
//*****************************************************************************
// Constructor
//*****************************************************************************
  template<class T>
TURBULENCE<T>::TURBULENCE(PARAMETERS<T> *par, MPI_DRIVER<T> *md, 
    CURVILINEAR_GRID<T> *g, ARRAY_3D<T> *rhoin, ARRAY_3D<VECTOR_3D<T> > *uin) 
: parameters(par), mpi_driver(md), grid(g), rho(rhoin), u(uin)
{
  imin = grid->I_Min(); imax = grid->I_Max(); 
  jmin = grid->J_Min(); jmax = grid->J_Max();
  kmin = grid->K_Min(); kmax = grid->K_Max();
  halo = grid->Halo_Size();

  eddy_viscosity = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
  eddy_diffusivity = new ARRAY_3D<T>(*eddy_viscosity);
  
  ss = new ARRAY_1D<ARRAY_3D<T>* >(9);
  tt = new ARRAY_1D<ARRAY_3D<T>* >(9);
  zz = new ARRAY_1D<ARRAY_3D<T>* >(9);
  for(m = 1; m <= 9, m++){
    (*ss)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
    (*tt)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
    (*zz)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
  }

  rr = new ARRAY_1D<ARRAY_3D<T>* >(6);
  for(m = 1; m <= 6, m++){
    (*rr)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
  }
}
//*****************************************************************************
// Destructor
//*****************************************************************************
  template<class T>
TURBULENCE<T>::~TURBULENCE()
{
  delete eddy_viscosity; delete eddy_diffusivity; delete tau;
  delete ss; delete tt; delete zz; delete rr;
}
//*****************************************************************************
// Calculate nu_t and kappa_t from Zang 1993 LES model
//*****************************************************************************
template <class T>
void Calculate_Eddy_Viscosity()
{
  tiny = (T)1e-20;
  alpha2 = (T)4;

  // T:
  // xi_x (1) et_x (2) zt_x (3)
  // xi_y (4) et_y (5) zt_y (6)
  // xi_z (7) et_z (8) zt_z (9)
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        

}
//*****************************************************************************
#endif
