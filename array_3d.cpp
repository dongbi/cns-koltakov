// array_3d.cpp
// 3D array class that takes care of hallo regions
#include "array_3d.h"
#include "vector_3d.h"
#include <cstring>

//****************************************************************************
template <class T>
void ARRAY_3D<T>::Init_Array(const int imin, const int imax, const int jmin,
      const int jmax, const int kmin, const int kmax, const int num_halo_cells)
{
  i_min=imin; i_max=imax; j_min=jmin; j_max=jmax; k_min=kmin; k_max=kmax;
  halo_size=num_halo_cells;
  i_size_w_h = i_max - i_min + 1 + 2*halo_size;
  j_size_w_h = j_max - j_min + 1 + 2*halo_size;
  k_size_w_h = k_max - k_min + 1 + 2*halo_size;
  // w/ halo regions
  i_min_w_h = -halo_size+i_min; i_max_w_h = i_max + halo_size;
  j_min_w_h = -halo_size+j_min; j_max_w_h = j_max + halo_size;
  k_min_w_h = -halo_size+k_min; k_max_w_h = k_max + halo_size;

  total_size = i_size_w_h*j_size_w_h*k_size_w_h;

  array = new T[total_size];
  memset(array, 0, total_size*sizeof(T));
}
//*****************************************************************************
template <class T>
void ARRAY_3D<T>::Set_All_Elements_To(const T& element_value)
{
  for(int n=0; n<total_size; n++) array[n] = element_value;
}
//*****************************************************************************
template <class T>
void ARRAY_3D<T>::Set_Elements_To(const T* a, int size)
{
  assert(size==total_size);
  for(int n=0; n<total_size; n++) array[n] = a[n];
}
//*****************************************************************************
template <class T>
void ARRAY_3D<T>::Set_NonHalo_Elements_To(const T& element_value)
{
  for(int i=I_Min();i<=I_Max();i++) 
    for(int j=J_Min();j<=J_Max();j++)  
      for(int k=K_Min();k<=K_Max();k++)
	(*this)(i,j,k) = element_value;
}
//*****************************************************************************
/*
template <class T>
T& ARRAY_3D<T>::operator () (const int i, const int j, const int k) 
{
  assert(i >= i_min_w_h && i <= i_max_w_h && j >= j_min_w_h && 
         j <= j_max_w_h && k >= k_min_w_h && k <= k_max_w_h);
  return array[i-i_min_w_h + (j-j_min_w_h)*i_size_w_h + 
		 (k-k_min_w_h)*i_size_w_h*j_size_w_h];
}

template <class T>
T ARRAY_3D<T>::operator () (const int i, const int j, const int k) const
{
  assert(i >= i_min_w_h && i <= i_max_w_h && j >= j_min_w_h && 
         j <= j_max_w_h && k >= k_min_w_h && k <= k_max_w_h);
  return array[i-i_min_w_h + (j-j_min_w_h)*i_size_w_h + 
		 (k-k_min_w_h)*i_size_w_h*j_size_w_h];
}
*/
template class ARRAY_3D<double>;
template class ARRAY_3D<VECTOR_3D<double> >;
