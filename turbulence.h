// turbulence.h

#ifndef __TURBULENCE__
#define __TURBULENCE__

#include "array_3d.h"

template <class T=double> 
class TURBULENCE
{
 public:
  TURBULENCE() : eddy_viscosity(NULL), eddy_diffusivity(NULL), tau(NULL) {}
  ~TURBULENCE() {}

  void Calculate_Eddy_Viscosity();

  ARRAY_3D<T> *eddy_viscosity, *eddy_diffusivity;
  ARRAY_3D<VECTOR_3D<T> > *tau;
};

#endif
