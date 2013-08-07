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

  void Filter(ARRAY_3D<T>& f, ARRAY_3D<T>& r, T weit, T fact);
  void SSBC(ARRAY_3D<T>& s);


  int i, j, k, m;
  int imin, imax, jmin, jmax, kmin, kmax, halo; //halo = 2
  int imin_w_h, imax_w_h, jmin_w_h, jmax_w_h, kmin_w_h, kmax_w_h; //halo = 1
  T tiny, alpha2, weit, fact, temp, fmin, fmax;
  ARRAY_1D<ARRAY_3D<T>* > *ss, *tt, *zz, *rr;
  ARRAY_3D<T> *sab;

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
  
  imin_w_h = grid->I_Min() + 1; imax_w_h = grid->I_Max() + 1; 
  jmin_w_h = grid->J_Min() + 1; jmax_w_h = grid->J_Max() + 1;
  kmin_w_h = grid->K_Min() + 1; kmax_w_h = grid->K_Max() + 1;

  eddy_viscosity = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
  eddy_diffusivity = new ARRAY_3D<T>(*eddy_viscosity);
  sab = new ARRAY_3D<T>(*eddy_viscosity);

  ss = new ARRAY_1D<ARRAY_3D<T>* >(9);
  tt = new ARRAY_1D<ARRAY_3D<T>* >(9);
  zz = new ARRAY_1D<ARRAY_3D<T>* >(9);
  for(m = 1; m <= 9; m++){
    (*ss)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
    (*tt)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
    (*zz)(m) = new ARRAY_3D<T>(imin,imax,jmin,jmax,kmin,kmax,halo-1);
  }

  rr = new ARRAY_1D<ARRAY_3D<T>* >(6);
  for(m = 1; m <= 6; m++){
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
  delete sab;
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
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        temp = (T)0.25 * grid->inverse_jacobian(i,j,k); 
        (*(*tt)(1))(i,j,k) = temp * ( grid->XI_X(i-1,j,k) + grid->XI_X(i,j,k) );
        (*(*tt)(2))(i,j,k) = temp * ( grid->XI_X(i,j-1,k) + grid->XI_X(i,j,k) );
        (*(*tt)(3))(i,j,k) = temp * ( grid->XI_X(i,j,k-1) + grid->XI_X(i,j,k) );
        (*(*tt)(4))(i,j,k) = temp * ( grid->XI_X(i-1,j,k) + grid->XI_X(i,j,k) );
        (*(*tt)(5))(i,j,k) = temp * ( grid->XI_X(i,j-1,k) + grid->XI_X(i,j,k) );
        (*(*tt)(6))(i,j,k) = temp * ( grid->XI_X(i,j,k-1) + grid->XI_X(i,j,k) );
        (*(*tt)(7))(i,j,k) = temp * ( grid->XI_X(i-1,j,k) + grid->XI_X(i,j,k) );
        (*(*tt)(8))(i,j,k) = temp * ( grid->XI_X(i,j-1,k) + grid->XI_X(i,j,k) );
        (*(*tt)(9))(i,j,k) = temp * ( grid->XI_X(i,j,k-1) + grid->XI_X(i,j,k) );

      }

  // Z:
  // uxi (1) uet (2) uzt (3)
  // vxi (4) vet (5) vzt (6)
  // wxi (7) wet (8) wzt (9)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*zz)(1))(i,j,k) = u(i+1,j,k).x - u(i-1,j,k).x;  
        (*(*zz)(2))(i,j,k) = u(i,j+1,k).x - u(i,j-1,k).x;  
        (*(*zz)(3))(i,j,k) = u(i,j,k+1).x - u(i,j,k-1).x;  
        (*(*zz)(4))(i,j,k) = u(i+1,j,k).y - u(i-1,j,k).y;  
        (*(*zz)(5))(i,j,k) = u(i,j+1,k).y - u(i,j-1,k).y;  
        (*(*zz)(6))(i,j,k) = u(i,j,k+1).y - u(i,j,k-1).y;  
        (*(*zz)(7))(i,j,k) = u(i+1,j,k).z - u(i-1,j,k).z;  
        (*(*zz)(8))(i,j,k) = u(i,j+1,k).z - u(i,j-1,k).z;  
        (*(*zz)(9))(i,j,k) = u(i,j,k+1).z - u(i,j,k-1).z;  
      }

  // R:
  // rxi (1) ret (2) rzt (3)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*rr)(1))(i,j,k) = rho(i+1,j,k) - rho(i-1,j,k);
        (*(*rr)(2))(i,j,k) = rho(i,j+1,k) - rho(i,j-1,k);
        (*(*rr)(3))(i,j,k) = rho(i,j,k+1) - rho(i,j,k-1);
      }

  // S:
  // S11 (1) S22 (2) S33 (3)
  // S12 (4) S23 (5) S31 (6)
  // rsx (7) rsy (8) rsz (9)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*ss)(1))(i,j,k) = (*(*tt)(1))(i,j,k) * (*(*zz)(1))(i,j,k)
                           + (*(*tt)(2))(i,j,k) * (*(*zz)(2))(i,j,k)
                           + (*(*tt)(3))(i,j,k) * (*(*zz)(3))(i,j,k);
        (*(*ss)(2))(i,j,k) = (*(*tt)(4))(i,j,k) * (*(*zz)(4))(i,j,k)
                           + (*(*tt)(5))(i,j,k) * (*(*zz)(5))(i,j,k)
                           + (*(*tt)(6))(i,j,k) * (*(*zz)(6))(i,j,k);
        (*(*ss)(3))(i,j,k) = (*(*tt)(7))(i,j,k) * (*(*zz)(7))(i,j,k)
                           + (*(*tt)(8))(i,j,k) * (*(*zz)(8))(i,j,k)
                           + (*(*tt)(9))(i,j,k) * (*(*zz)(9))(i,j,k);
        (*(*ss)(4))(i,j,k) = (*(*tt)(4))(i,j,k) * (*(*zz)(1))(i,j,k)
                           + (*(*tt)(5))(i,j,k) * (*(*zz)(2))(i,j,k)
                           + (*(*tt)(6))(i,j,k) * (*(*zz)(3))(i,j,k)
                           + (*(*tt)(1))(i,j,k) * (*(*zz)(4))(i,j,k)
                           + (*(*tt)(2))(i,j,k) * (*(*zz)(5))(i,j,k)
                           + (*(*tt)(3))(i,j,k) * (*(*zz)(6))(i,j,k);
        (*(*ss)(5))(i,j,k) = (*(*tt)(7))(i,j,k) * (*(*zz)(4))(i,j,k)
                           + (*(*tt)(8))(i,j,k) * (*(*zz)(5))(i,j,k)
                           + (*(*tt)(9))(i,j,k) * (*(*zz)(6))(i,j,k)
                           + (*(*tt)(4))(i,j,k) * (*(*zz)(7))(i,j,k)
                           + (*(*tt)(5))(i,j,k) * (*(*zz)(8))(i,j,k)
                           + (*(*tt)(6))(i,j,k) * (*(*zz)(9))(i,j,k);
        (*(*ss)(6))(i,j,k) = (*(*tt)(1))(i,j,k) * (*(*zz)(7))(i,j,k)
                           + (*(*tt)(2))(i,j,k) * (*(*zz)(8))(i,j,k)
                           + (*(*tt)(3))(i,j,k) * (*(*zz)(9))(i,j,k)
                           + (*(*tt)(7))(i,j,k) * (*(*zz)(1))(i,j,k)
                           + (*(*tt)(8))(i,j,k) * (*(*zz)(2))(i,j,k)
                           + (*(*tt)(9))(i,j,k) * (*(*zz)(3))(i,j,k);
        (*(*ss)(4))(i,j,k) *= (T)0.5;
        (*(*ss)(5))(i,j,k) *= (T)0.5;
        (*(*ss)(6))(i,j,k) *= (T)0.5;
        (*(*ss)(7))(i,j,k) = (*(*tt)(1))(i,j,k) * (*(*rr)(1))(i,j,k)
                           + (*(*tt)(2))(i,j,k) * (*(*rr)(2))(i,j,k)
                           + (*(*tt)(3))(i,j,k) * (*(*rr)(3))(i,j,k);
        (*(*ss)(7))(i,j,k) = (*(*tt)(4))(i,j,k) * (*(*rr)(1))(i,j,k)
                           + (*(*tt)(5))(i,j,k) * (*(*rr)(2))(i,j,k)
                           + (*(*tt)(6))(i,j,k) * (*(*rr)(3))(i,j,k);
        (*(*ss)(7))(i,j,k) = (*(*tt)(7))(i,j,k) * (*(*rr)(1))(i,j,k)
                           + (*(*tt)(8))(i,j,k) * (*(*rr)(2))(i,j,k)
                           + (*(*tt)(9))(i,j,k) * (*(*rr)(3))(i,j,k);

      }
  // |S|
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        sab(i,j,k) = pow((*(*ss)(1))(i,j,k),2) 
                   + pow((*(*ss)(2))(i,j,k),2) 
                   + pow((*(*ss)(3))(i,j,k),2) 
             + 2 * ( pow((*(*ss)(4))(i,j,k),2) 
                   + pow((*(*ss)(5))(i,j,k),2) 
                   + pow((*(*ss)(6))(i,j,k),2) );
        sab(i,j,k) = sqrt(2*sab(i,j,k));
      }

  // Z:
  // |S|S11 (1) |S|S22 (2) |S|S33 (3)
  // |S|S12 (4) |S|S23 (5) |S|S31 (6)
  // |S|rsx (7) |S|rsy (8) |S|rsz (9)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*zz)(1))(i,j,k) = sab(i,j,k) * (*(*ss)(1))(i,j,k);
        (*(*zz)(2))(i,j,k) = sab(i,j,k) * (*(*ss)(2))(i,j,k);
        (*(*zz)(3))(i,j,k) = sab(i,j,k) * (*(*ss)(3))(i,j,k);
        (*(*zz)(4))(i,j,k) = sab(i,j,k) * (*(*ss)(4))(i,j,k);
        (*(*zz)(5))(i,j,k) = sab(i,j,k) * (*(*ss)(5))(i,j,k);
        (*(*zz)(6))(i,j,k) = sab(i,j,k) * (*(*ss)(6))(i,j,k);
        (*(*zz)(7))(i,j,k) = sab(i,j,k) * (*(*ss)(7))(i,j,k);
        (*(*zz)(8))(i,j,k) = sab(i,j,k) * (*(*ss)(8))(i,j,k);
        (*(*zz)(9))(i,j,k) = sab(i,j,k) * (*(*ss)(9))(i,j,k);
      }
  
  // Z:
  // <|S|S11> (1) <|S|S22> (2) <|S|S33> (3)
  // <|S|S12> (4) <|S|S23> (5) <|S|S31> (6)
  // <|S|rsx> (7) <|S|rsy> (8) <|S|rsz> (9)
  weit = (T)2;
  fact = (T)1/64;
  Filter(*(*zz)(1), *(*rr)(5), weit, fact);
  Filter(*(*zz)(2), *(*rr)(5), weit, fact);
  Filter(*(*zz)(3), *(*rr)(5), weit, fact);
  Filter(*(*zz)(4), *(*rr)(5), weit, fact);
  Filter(*(*zz)(5), *(*rr)(5), weit, fact);
  Filter(*(*zz)(6), *(*rr)(5), weit, fact);
  Filter(*(*zz)(7), *(*rr)(5), weit, fact);
  Filter(*(*zz)(8), *(*rr)(5), weit, fact);
  Filter(*(*zz)(9), *(*rr)(5), weit, fact);

  // R:
  // |S| (4)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*rr)(4))(i,j,k) = sab(i,j,k);
      }

  // R:
  // |S| (4)
  weit = (T)2;
  fact = (T)1/64;
  Filter((*rr)(4), (*rr)(5), weit, fact);

  // S:
  // <S11> (1) <S22> (2) <S33> (3)
  // <S12> (4) <S23> (5) <S31> (6)
  // <rsx> (7) <rsy> (8) <rsz> (9)
  weit = (T)2;
  fact = (T)1/64;
  Filter(*(*ss)(1), *(*rr)(5), weit, fact);
  Filter(*(*ss)(2), *(*rr)(5), weit, fact);
  Filter(*(*ss)(3), *(*rr)(5), weit, fact);
  Filter(*(*ss)(4), *(*rr)(5), weit, fact);
  Filter(*(*ss)(5), *(*rr)(5), weit, fact);
  Filter(*(*ss)(6), *(*rr)(5), weit, fact);
  Filter(*(*ss)(7), *(*rr)(5), weit, fact);
  Filter(*(*ss)(8), *(*rr)(5), weit, fact);
  Filter(*(*ss)(9), *(*rr)(5), weit, fact);
  
  // Z:
  // M11 (1) M22 (2) M33 (3)
  // M12 (4) M23 (5) M31 (6)
  // N1  (7) N2  (8) N3  (9)
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*zz)(1))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(1)(i,j,k) - (*(*zz)(1))(i,j,k);
        (*(*zz)(2))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(2)(i,j,k) - (*(*zz)(2))(i,j,k);
        (*(*zz)(3))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(3)(i,j,k) - (*(*zz)(3))(i,j,k);
        (*(*zz)(4))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(4)(i,j,k) - (*(*zz)(4))(i,j,k);
        (*(*zz)(5))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(5)(i,j,k) - (*(*zz)(5))(i,j,k);
        (*(*zz)(6))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(6)(i,j,k) - (*(*zz)(6))(i,j,k);
        (*(*zz)(7))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(7)(i,j,k) - (*(*zz)(7))(i,j,k);
        (*(*zz)(8))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(8)(i,j,k) - (*(*zz)(8))(i,j,k);
        (*(*zz)(9))(i,j,k) = alpha2 * (*(*rr)(4))(i,j,k) * (*(*ss)(9)(i,j,k) - (*(*zz)(9))(i,j,k);
      }

  // T:
  // uu (1) vv (2) ww (3)
  // uv (4) vw (5) wu (6)
  // ru (7) rv (8) rw (9)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*(tt)(1))(i,j,k) = pow(u(i,j,k).x,2);
        (*(*(tt)(2))(i,j,k) = pow(u(i,j,k).y,2);
        (*(*(tt)(3))(i,j,k) = pow(u(i,j,k).z,2);
        (*(*(tt)(4))(i,j,k) = u(i,j,k).x * u(i,j,k).y; 
        (*(*(tt)(5))(i,j,k) = u(i,j,k).y * u(i,j,k).z;
        (*(*(tt)(6))(i,j,k) = u(i,j,k).z * u(i,j,k).x;
        (*(*(tt)(7))(i,j,k) = phi(i,j,k) * u(i,j,k).x;
        (*(*(tt)(8))(i,j,k) = phi(i,j,k) * u(i,j,k).y;
        (*(*(tt)(9))(i,j,k) = phi(i,j,k) * u(i,j,k).z;
      }

  // T:
  // <uu> (1) <vv> (2) <ww> (3)
  // <uv> (4) <vw> (5) <wu> (6)
  // <ru> (7) <rv> (8) <rw> (9)
  weit = (T)2;
  fact = (T)1/64;
  Filter(*(*tt)(1), *(*rr)(5), weit, fact);
  Filter(*(*tt)(2), *(*rr)(5), weit, fact);
  Filter(*(*tt)(3), *(*rr)(5), weit, fact);
  Filter(*(*tt)(4), *(*rr)(5), weit, fact);
  Filter(*(*tt)(5), *(*rr)(5), weit, fact);
  Filter(*(*tt)(6), *(*rr)(5), weit, fact);
  Filter(*(*tt)(7), *(*rr)(5), weit, fact);
  Filter(*(*tt)(8), *(*rr)(5), weit, fact);
  Filter(*(*tt)(9), *(*rr)(5), weit, fact);

  // R:
  // u (1) v (2) w (3) r (4)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*rr)(1))(i,j,k) = u(i,j,k).x;
        (*(*rr)(2))(i,j,k) = u(i,j,k).y;
        (*(*rr)(3))(i,j,k) = u(i,j,k).z;
        (*(*rr)(4))(i,j,k) = phi(i,j,k);
      }

  // R:
  // <u> (1) <v> (2) <w> (3) <r> (4);
  weit = (T)2;
  fact = (T)1/64;
  Filter(*(*rr)(1), *(*rr)(5), weit, fact);
  Filter(*(*rr)(2), *(*rr)(5), weit, fact);
  Filter(*(*rr)(3), *(*rr)(5), weit, fact);
  Filter(*(*rr)(4), *(*rr)(5), weit, fact);

  // T:
  // L11 (1) L22 (2) L33 (3)
  // L12 (4) L23 (5) L31 (6)
  // K1  (7) K2  (8) K3  (9)
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*tt)(1))(i,j,k) -= pow((*(*rr)(1))(i,j,k),2); 
        (*(*tt)(2))(i,j,k) -= pow((*(*rr)(2))(i,j,k),2); 
        (*(*tt)(3))(i,j,k) -= pow((*(*rr)(3))(i,j,k),2); 
        (*(*tt)(4))(i,j,k) -= (*(*rr)(1))(i,j,k) * (*(*rr)(2))(i,j,k);
        (*(*tt)(5))(i,j,k) -= (*(*rr)(2))(i,j,k) * (*(*rr)(3))(i,j,k);
        (*(*tt)(6))(i,j,k) -= (*(*rr)(3))(i,j,k) * (*(*rr)(1))(i,j,k);
        (*(*tt)(7))(i,j,k) -= (*(*rr)(4))(i,j,k) * (*(*rr)(1))(i,j,k);
        (*(*tt)(8))(i,j,k) -= (*(*rr)(4))(i,j,k) * (*(*rr)(2))(i,j,k);
        (*(*tt)(9))(i,j,k) -= (*(*rr)(4))(i,j,k) * (*(*rr)(3))(i,j,k);
      }
  
  // R:
  // u (1) v (2) w (3) r (4)
  for(i = imin_w_h; i <= imax_w_h; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        (*(*rr)(1))(i,j,k) = u(i,j,k).x;
        (*(*rr)(2))(i,j,k) = u(i,j,k).y;
        (*(*rr)(3))(i,j,k) = u(i,j,k).z;
        (*(*rr)(4))(i,j,k) = phi(i,j,k);
      }
  
  // R:
  // [u] (1) [v] (2) [w] (3) [r] (4);
  weit = (T)6;
  fact = (T)1/512;
  Filter(*(*rr)(1), *(*rr)(5), weit, fact);
  Filter(*(*rr)(2), *(*rr)(5), weit, fact);
  Filter(*(*rr)(3), *(*rr)(5), weit, fact);
  Filter(*(*rr)(4), *(*rr)(5), weit, fact);
  SSBC(*(*rr)(1));
  SSBC(*(*rr)(2));
  SSBC(*(*rr)(3));
  SSBC(*(*rr)(4));
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*(*rr)(1),1);
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*(*rr)(2),1);
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*(*rr)(3),1);
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*(*rr)(4),1);

  // S:
  // [u][u] (1) [v][v] (2) [w][w] (3)
  // [u][v] (4) [v][w] (5) [w][u] (6)
  // [r][u] (7) [r][v] (8) [r][w] (9)
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*ss)(1))(i,j,k) = pow((*(*rr)(1))(i,j,k),2); 
        (*(*ss)(2))(i,j,k) = pow((*(*rr)(2))(i,j,k),2); 
        (*(*ss)(3))(i,j,k) = pow((*(*rr)(3))(i,j,k),2); 
        (*(*ss)(4))(i,j,k) = (*(*rr)(1))(i,j,k) * (*(*rr)(2))(i,j,k);
        (*(*ss)(5))(i,j,k) = (*(*rr)(2))(i,j,k) * (*(*rr)(3))(i,j,k);
        (*(*ss)(6))(i,j,k) = (*(*rr)(3))(i,j,k) * (*(*rr)(1))(i,j,k);
        (*(*ss)(7))(i,j,k) = (*(*rr)(4))(i,j,k) * (*(*rr)(1))(i,j,k);
        (*(*ss)(8))(i,j,k) = (*(*rr)(4))(i,j,k) * (*(*rr)(2))(i,j,k);
        (*(*ss)(9))(i,j,k) = (*(*rr)(4))(i,j,k) * (*(*rr)(3))(i,j,k);
      }

  // S:
  // <[u][u]> (1) <[v][v]> (2) <[w][w]> (3)
  // <[u][v]> (4) <[v][w]> (5) <[w][u]> (6)
  // <[r][u]> (7) <[r][v]> (8) <[r][w]> (9)
  weit = (T)2;
  fact = (T)1/64;
  Filter(*(*ss)(1), *(*rr)(5), weit, fact);
  Filter(*(*ss)(2), *(*rr)(5), weit, fact);
  Filter(*(*ss)(3), *(*rr)(5), weit, fact);
  Filter(*(*ss)(4), *(*rr)(5), weit, fact);
  Filter(*(*ss)(5), *(*rr)(5), weit, fact);
  Filter(*(*ss)(6), *(*rr)(5), weit, fact);
  Filter(*(*ss)(7), *(*rr)(5), weit, fact);
  Filter(*(*ss)(8), *(*rr)(5), weit, fact);
  Filter(*(*ss)(9), *(*rr)(5), weit, fact);

  // R:
  // <[u]> (1) <[v]> (2) <[w]> (3) <[r]> (4);
  weit = (T)2;
  fact = (T)1/64;
  Filter(*(*rr)(1), *(*rr)(5), weit, fact);
  Filter(*(*rr)(2), *(*rr)(5), weit, fact);
  Filter(*(*rr)(3), *(*rr)(5), weit, fact);
  Filter(*(*rr)(4), *(*rr)(5), weit, fact);

  // S:
  // H11 (1) H22 (2) H33 (3)
  // H12 (4) H23 (5) H31 (6)
  // J1  (7) J2  (8) J3  (9)
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*ss)(1))(i,j,k) -= pow((*(*rr)(1))(i,j,k),2); 
        (*(*ss)(2))(i,j,k) -= pow((*(*rr)(2))(i,j,k),2); 
        (*(*ss)(3))(i,j,k) -= pow((*(*rr)(3))(i,j,k),2); 
        (*(*ss)(4))(i,j,k) -= (*(*rr)(1))(i,j,k) * (*(*rr)(2))(i,j,k);
        (*(*ss)(5))(i,j,k) -= (*(*rr)(2))(i,j,k) * (*(*rr)(3))(i,j,k);
        (*(*ss)(6))(i,j,k) -= (*(*rr)(3))(i,j,k) * (*(*rr)(1))(i,j,k);
        (*(*ss)(7))(i,j,k) -= (*(*rr)(4))(i,j,k) * (*(*rr)(1))(i,j,k);
        (*(*ss)(8))(i,j,k) -= (*(*rr)(4))(i,j,k) * (*(*rr)(2))(i,j,k);
        (*(*ss)(9))(i,j,k) -= (*(*rr)(4))(i,j,k) * (*(*rr)(3))(i,j,k);
      }

  // MijMij => R1
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(1))(i,j,k) = pow((*(*zz)(1))(i,j,k),2) 
                           + pow((*(*zz)(2))(i,j,k),2) 
                           + pow((*(*zz)(3))(i,j,k),2) 
                     + 2 * ( pow((*(*zz)(4))(i,j,k),2) 
                           + pow((*(*zz)(5))(i,j,k),2) 
                           + pow((*(*zz)(6))(i,j,k),2) );
      }

  // MijLij => R2
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(2))(i,j,k) = (*(*zz)(1))(i,j,k) * (*(*tt)(1))(i,j,k)
                           + (*(*zz)(2))(i,j,k) * (*(*tt)(2))(i,j,k)
                           + (*(*zz)(3))(i,j,k) * (*(*tt)(3))(i,j,k)
                     + 2 * ( (*(*zz)(4))(i,j,k) * (*(*tt)(4))(i,j,k)
                           + (*(*zz)(5))(i,j,k) * (*(*tt)(5))(i,j,k)
                           + (*(*zz)(6))(i,j,k) * (*(*tt)(6))(i,j,k) );
      }

  // MijHij => R3
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(3))(i,j,k) = (*(*zz)(1))(i,j,k) * (*(*ss)(1))(i,j,k)
                           + (*(*zz)(2))(i,j,k) * (*(*ss)(2))(i,j,k)
                           + (*(*zz)(3))(i,j,k) * (*(*ss)(3))(i,j,k)
                     + 2 * ( (*(*zz)(4))(i,j,k) * (*(*ss)(4))(i,j,k)
                           + (*(*zz)(5))(i,j,k) * (*(*ss)(5))(i,j,k)
                           + (*(*zz)(6))(i,j,k) * (*(*ss)(6))(i,j,k) );
      }

  // MijMij => Z1
  // MijLij => Z2
  // MijHij => Z3
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*zz)(1))(i,j,k) = (*(*rr)(1))(i,j,k);
        (*(*zz)(2))(i,j,k) = (*(*rr)(2))(i,j,k);
        (*(*zz)(3))(i,j,k) = (*(*rr)(3))(i,j,k);
      }

  //NiNi => R1
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(1))(i,j,k) = pow((*(*zz)(7))(i,j,k),2)
                           + pow((*(*zz)(8))(i,j,k),2)
                           + pow((*(*zz)(9))(i,j,k),2);
      }

  //NiKi => R2
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(2))(i,j,k) = (*(*zz)(7))(i,j,k) * (*(*tt)(7))(i,j,k)
                           + (*(*zz)(8))(i,j,k) * (*(*tt)(8))(i,j,k)
                           + (*(*zz)(9))(i,j,k) * (*(*tt)(9))(i,j,k);
      }

  //NiJi => R3
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(3))(i,j,k) = (*(*zz)(7))(i,j,k) * (*(*ss)(7))(i,j,k)
                           + (*(*zz)(8))(i,j,k) * (*(*ss)(8))(i,j,k)
                           + (*(*zz)(9))(i,j,k) * (*(*ss)(9))(i,j,k);
      }

  //NiNi => Z4
  //NiKi => Z5
  //NiJi => Z6
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*zz)(4))(i,j,k) = (*(*rr)(1))(i,j,k);
        (*(*zz)(5))(i,j,k) = (*(*rr)(2))(i,j,k);
        (*(*zz)(6))(i,j,k) = (*(*rr)(3))(i,j,k);
      }

  //               Z2       Z3
  //           - MijLij + MijHij
  // 2*C*D^2 = -------------------
  //                  MijMij
  //                    Z1

  //                Z5     Z6
  //             - NiKi + NiJi
  // C*D^2/Prt = ---------------
  //                  NiNi
  //                   Z4

  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*zz)(1))(i,j,k) = max((*(*zz)(1))(i,j,k), tiny);
        (*(*zz)(4))(i,j,k) = max((*(*zz)(4))(i,j,k), tiny);
      }

  // C*D^2 => R1
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(1))(i,j,k) = -(*(*zz)(2))(i,j,k) + (*(*zz)(3))(i,j,k);
        (*(*rr)(1))(i,j,k) *= (T)0.5 / (*(*zz)(1))(i,j,k);
      }

  // C/Prt*D^2 => R2
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(2))(i,j,k) = -(*(*zz)(5))(i,j,k) + (*(*zz)(6))(i,j,k);
        (*(*rr)(2))(i,j,k) /= (*(*zz)(4))(i,j,k);
      }

  // 1/D^2 => R3
  temp = (T)2/3;
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(3))(i,j,k) = pow(grid->inverse_jacobian(i,j,k),temp);
      }

  // C => R2
  // C/Prt => R3
  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        (*(*rr)(1))(i,j,k) *= (*(*rr)(3))(i,j,k);
        (*(*rr)(2))(i,j,k) *= (*(*rr)(3))(i,j,k);
      }
  SSBC(*(*rr)(1));
  SSBC(*(*rr)(2));
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*(*rr)(1),1);
  mpi_driver->Exchange_Ghost_Values_For_Scalar_Field(*(*rr)(2),1);

}
//*****************************************************************************
// Spatial filter
//*****************************************************************************
  template<class T>
void Filter(ARRAY_3D<T>& f, ARRAY_3D<T>& r, T weit, T fact)
{
  int i, j, k;

  for(i = imin; i <= imax; i++)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        r(i,j,k) = f(i-1,j,k) + weit * f(i,j,k) + f(i+1,j,k);
      }

  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        f(i,j,k) = r(i,j-1,j) + weit * r(i,j,k) + r(i,j+1,k);
      }

  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        r(i,j,k) = f(i,j,k-1) + weit * f(i,j,k) + f(i,j,k+1);
      }

  for(i = imin; i <= imax; i++)
    for(j = jmin; j <= jmax; j++)
      for(k = kmin; k <= kmax; k++){
        f(i,j,k) = fact * r(i,j,k); 
      }

}
//*****************************************************************************
// Boundary Condition
//*****************************************************************************
  template<class T>
void SSBC(ARRAY_3D<T>& s)
{
  int i, j, k;

  //West BC
  if(mpi_driver->west_proc == MPI_PROC_NULL)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        s(imin_w_h,j,k) = (T)2 * s(imin,j,k) - s(imin+1,j,k);
      }

  //East BC
  if(mpi_driver->east_proc == MPI_PROC_NULL)
    for(j = jmin_w_h; j <= jmax_w_h; j++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        s(imax_w_h,j,k) = (T)2 * s(imax,j,k) - s(imax-1,j,k);
      }

  //South BC
  if(mpi_driver->suth_proc == MPI_PROC_NULL)
    for(i = imin_w_h; i <= imax_w_h; i++)
      for(j = jmin_w_h; j <= jmax_w_h; j++){
        s(i,j,kmin_w_h) = (T)2 * s(i,j,kmin) - s(imin+1,j,k);
      }

  //North BC
  if(mpi_driver->nrth_proc == MPI_PROC_NULL)
    for(i = imin_w_h; i <= imax_w_h; i++)
      for(j = jmin_w_h; j <= jmax_w_h; j++){
        s(i,j,kmax_w_h) = (T)2 * s(i,j,kmax) - s(i,j,kmax-1);
      }

  //Front BC
  if(mpi_driver->frnt_proc == MPI_PROC_NULL)
    for(i = imin_w_h; i <= imax_w_h; i++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        s(i,jmin_w_h,k) = (T)2 * s(i,jmin,k) - s(i,jmin+1,k);
      }

  //Back BC
  if(mpi_driver->back_proc == MPI_PROC_NULL)
    for(i = imin_w_h; i <= imax_w_h; i++)
      for(k = kmin_w_h; k <= kmax_w_h; k++){
        s(i,jmax_w_h,k) = (T)2 * s(i,jmax,k) - s(i,jmax-1,k);
      }
}
//*****************************************************************************
#endif
