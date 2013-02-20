// array_2d.h
//
// 2D array class that takes care of hallo regions
// Index(for x,y): -halo_size+i_min...i_min...i_max...i_max+halo_size
// all non-ghost indexes start with i_min and end with i_max
//
#ifndef __ARRAY_2D__
#define __ARRAY_2D__
#include <iostream>
#include <assert.h>
#include <cstring>
//#include "array_3d.h"

using namespace std;

template <class T> class ARRAY_3D;

template <class T=double> 
class ARRAY_2D
{
 public:
  ARRAY_2D(const int imin, const int imax, const int jmin, const int jmax,
	   const int num_halo_cells)
  {Init_Array(imin, imax, jmin, jmax, num_halo_cells);}

  ARRAY_2D(const ARRAY_2D<T>& a)
  {Init_Array(a.I_Min(), a.I_Max(), a.J_Min(), a.J_Max(), a.Halo_Size());
   for(int n=0;n<total_size;n++) array[n] = a.Raw_Array(n);}
  
  template<class T2>
  ARRAY_2D(const ARRAY_2D<T2>& at)
  {Init_Array(at.I_Min(), at.I_Max(), at.J_Min(), at.J_Max(), at.Halo_Size());}

  ARRAY_2D(const ARRAY_3D<T>& a3d, const int drop_component, 
                                   const int index_in_dropped_dimension = 1);
  ~ARRAY_2D() { delete [] array;}

  void Init_Array(const int imin, const int imax, 
		  const int jmin, const int jmax, const int num_halo_cells);
  void Set_All_Elements_To(const T element_value);

  T& operator() (const int i, const int j);
  T  operator() (const int i, const int j) const;
 
  ARRAY_2D<T>& operator=(const ARRAY_2D<T>& source)
  {Equal_Dimensions(*this,source);
   for(int n=0;n<total_size;n++) array[n] = source.Raw_Array(n); return *this;}

  ARRAY_2D<T>& operator+=(const ARRAY_2D<T>& source)
  {Equal_Dimensions(*this,source);
   for(int n=0;n<total_size;n++) array[n] += source.Raw_Array(n); return *this;}

  ARRAY_2D<T>& operator-=(const ARRAY_2D<T>& source)
  {Equal_Dimensions(*this,source);
   for(int n=0;n<total_size;n++) array[n] -= source.Raw_Array(n); return *this;}

  ARRAY_2D<T>& operator*=(const ARRAY_2D<T>& source) //element-wise
  {Equal_Dimensions(*this,source);
   for(int n=0;n<total_size;n++) array[n] *= source.Raw_Array(n); return *this;}

  ARRAY_2D<T>& operator/=(const ARRAY_2D<T>& source) //element-wise
  {Equal_Dimensions(*this,source);
   for(int n=0;n<total_size;n++) array[n] /= source.Raw_Array(n); return *this;}

  ARRAY_2D<T> operator*(const T scalar)
  {ARRAY_2D<T> a_copy(*this);
   for(int n=0;n<total_size;n++) a_copy.Raw_Array(n) *= scalar; return a_copy;}

  template<class T2>
  static bool Equal_Dimensions(const ARRAY_2D<T>& a, const ARRAY_2D<T2>& b)
  {return a.I_Min()==b.I_Min() && a.I_Max()==b.I_Max() && a.J_Min()==b.J_Min()
       && a.J_Max()==b.J_Max() && a.Halo_Size()==b.Halo_Size();}

  int I_Min() const {return i_min;}
  int I_Max() const {return i_max;}
  int J_Min() const {return j_min;}
  int J_Max() const {return j_max;}

  int I_Min_With_Halo() const {return i_min_w_h;}
  int I_Max_With_Halo() const {return i_max_w_h;}
  int J_Min_With_Halo() const {return j_min_w_h;}
  int J_Max_With_Halo() const {return j_max_w_h;}

  int I_Size_With_Halo() const {return i_size_w_h;}
  int J_Size_With_Halo() const {return j_size_w_h;}
  int Halo_Size() const {return halo_size;}

  T& Raw_Array(int n) {return array[n];}
  T Raw_Array(int n) const {return array[n];}

 private:
  int i_min, i_max, j_min, j_max; // w/o halo regions
  int i_min_w_h, i_max_w_h, j_min_w_h, j_max_w_h;
  int i_size_w_h, j_size_w_h, halo_size;
  int total_size;
  T* array;
};
//***************************************************************************
// Create 2D array from 3D array by dropping 1 component:
// x->drop_component = 1, y->drop_component = 2, z->drop_component = 3
//***************************************************************************
template <class T>
ARRAY_2D<T>::ARRAY_2D(const ARRAY_3D<T>& a3d, const int drop_component, 
                      const int index_in_dropped_dimension)
{
  assert(drop_component >= 1 && drop_component <= 3);
  int imin, imax, jmin, jmax;
  switch(drop_component){
   case 1: //use y,z
    imin=a3d.J_Min();imax=a3d.J_Max();jmin=a3d.K_Min();jmax=a3d.K_Max();break;
   case 2: //uze x,z
    imin=a3d.I_Min();imax=a3d.I_Max();jmin=a3d.K_Min();jmax=a3d.K_Max();break;
   case 3: //uze z,y
    imin=a3d.I_Min();imax=a3d.I_Max();jmin=a3d.J_Min();jmax=a3d.J_Max();break;
  }
  Init_Array(imin, imax, jmin, jmax, a3d.Halo_Size());
  //copy part of source array
  for(int i=imin; i<=imax; i++) 
    for(int j=jmin; j<=jmax; j++){
      T element;
      switch(drop_component){
       case 1: element = a3d(index_in_dropped_dimension, i, j);break;
       case 2: element = a3d(i, index_in_dropped_dimension, j);break;
       case 3: element = a3d(i, j, index_in_dropped_dimension);break;
      }				  
      (*this)(i,j) = element;
    }
}
//***************************************************************************
template <class T>
void ARRAY_2D<T>::Init_Array(const int imin, const int imax, const int jmin, 
                             const int jmax, const int num_halo_cells)
{
  i_min=imin; i_max=imax; j_min=jmin; j_max=jmax; halo_size=num_halo_cells;
  i_size_w_h = i_max - i_min + 1 + 2*halo_size;
  j_size_w_h = j_max - j_min + 1 + 2*halo_size;
  // w/ halo regions
  i_min_w_h = -halo_size+i_min; i_max_w_h = i_max + halo_size;
  j_min_w_h = -halo_size+j_min; j_max_w_h = j_max + halo_size;
  total_size = i_size_w_h * j_size_w_h;

  array = new T[total_size];
  //Set_All_Elements_To((T)0);
  memset(array, 0, total_size*sizeof(T));
}
//***************************************************************************
template <class T>
void ARRAY_2D<T>::Set_All_Elements_To(const T element_value)
{
  for(int n=0; n<total_size; n++) array[n] = element_value;
}
//***************************************************************************
template <class T>
T& ARRAY_2D<T>::operator () (const int i, const int j) 
{
  assert(i >= i_min_w_h && i <= i_max_w_h && j >= j_min_w_h && j <= j_max_w_h);
  return array[i-i_min_w_h + (j-j_min_w_h)*i_size_w_h];
}
//***************************************************************************
template <class T>
T ARRAY_2D<T>::operator () (const int i, const int j) const
{
  assert(i >= i_min_w_h && i <= i_max_w_h && j >= j_min_w_h && j <= j_max_w_h);
  return array[i-i_min_w_h + (j-j_min_w_h)*i_size_w_h];
}
//***************************************************************************
// global function
template <class T>
inline ostream& operator<< (ostream& output, const ARRAY_2D<T>& array)
{
  for (int i = array.I_Min_With_Halo(); i <= array.I_Max_With_Halo(); i++)
    for (int j = array.J_Min_With_Halo(); j <= array.J_Max_With_Halo(); j++)
	output <<"["<<i<<","<<j<<"]"<<"="<<array(i,j)<<endl;
  return output;
}
//***************************************************************************
#endif
