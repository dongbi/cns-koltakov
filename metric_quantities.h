//*****************************************************************************
// Metric Quantities class
//*****************************************************************************
#ifndef __METRIC_QUANTITIES__
#define __METRIC_QUANTITIES__

#include "array_3d.h"

template <class T=double>
class METRIC_QUANTITIES
{
 public: 
  METRIC_QUANTITIES() : G11(0),G12(0),G13(0),G21(0),G22(0),G23(0),
                        G31(0),G32(0),G33(0),GCC(0),delete_arrays(false)
  {}
  ~METRIC_QUANTITIES() 
  {
    if(delete_arrays) {
      delete G11; delete G12; delete G13; delete G21; delete G22; delete G23; 
      delete G31; delete G32; delete G33; delete GCC;
    }
  }
  void Set_Pointers(ARRAY_3D<T> *g11, ARRAY_3D<T> *g12, ARRAY_3D<T> *g13, 
                    ARRAY_3D<T> *g21, ARRAY_3D<T> *g22, ARRAY_3D<T> *g23,
         ARRAY_3D<T> *g31, ARRAY_3D<T> *g32, ARRAY_3D<T> *g33, ARRAY_3D<T> *gcc)
  {
    delete_arrays=false; G11=g11; G12=g12; G13=g13; G21=g21; G22=g22; G23=g23;
    G31=g31; G32=g32; G33=g33; GCC=gcc;
  }
  void Create_Copies(ARRAY_3D<T> *g11, ARRAY_3D<T> *g12, ARRAY_3D<T> *g13, 
                    ARRAY_3D<T> *g21, ARRAY_3D<T> *g22, ARRAY_3D<T> *g23,
	 ARRAY_3D<T> *g31, ARRAY_3D<T> *g32, ARRAY_3D<T> *g33, ARRAY_3D<T> *gcc)
  {
    delete_arrays=true;
    G11 = new ARRAY_3D<T>(*g11); G12 = new ARRAY_3D<T>(*g12); 
    G13 = new ARRAY_3D<T>(*g13); G21 = new ARRAY_3D<T>(*g21);
    G22 = new ARRAY_3D<T>(*g22); G23 = new ARRAY_3D<T>(*g23); 
    G31 = new ARRAY_3D<T>(*g31); G32 = new ARRAY_3D<T>(*g32);
    G33 = new ARRAY_3D<T>(*g33); GCC = new ARRAY_3D<T>(*gcc);
  }

  ARRAY_3D<T> *G11, *G12, *G13, *G21, *G22, *G23, *G31, *G32, *G33, *GCC;

 private:
  bool delete_arrays;
};
//*****************************************************************************
#endif
