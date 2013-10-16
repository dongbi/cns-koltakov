// array_3d.h
//
// 3D array class that takes care of hallo regions
// Index(for x,y,z): -halo+i_min...i_min...i_max-halo...i_max
// all non-ghost indexes start with i/j/k_min and end with i/j/k_max
//
#ifndef __ARRAY_3D__
#define __ARRAY_3D__
#include <iostream>
#include <assert.h>
#include <cmath>
#include "array_2d.h"
#include "vector_3d.h"

using namespace std;

template <class T=double> 
class ARRAY_3D
{
 public:
  ARRAY_3D(const int imin, const int imax, const int jmin, const int jmax,
	   const int kmin, const int kmax, const int num_halo_cells)
  {Init_Array(imin, imax, jmin, jmax, kmin, kmax, num_halo_cells);}
  
  ARRAY_3D(const ARRAY_3D<T>& ia, bool copy_contents = true)
  {Init_Array(ia.I_Min(), ia.I_Max(), ia.J_Min(), ia.J_Max(), 
		    ia.K_Min(), ia.K_Max(), ia.Halo_Size());
   if(copy_contents) for(int n=0;n<total_size;n++) array[n] = ia.Raw_Array(n);}
  
  template <class T2>
  ARRAY_3D(const ARRAY_3D<T2>& ta)
  {Init_Array(ta.I_Min(), ta.I_Max(), ta.J_Min(), ta.J_Max(), 
	      ta.K_Min(), ta.K_Max(), ta.Halo_Size());}

  ARRAY_3D(const ARRAY_2D<T>& i2da, int k_size=1) // 3D array w/ same halo
  {Init_Array(i2da.I_Min(), i2da.I_Max(), 
	      i2da.J_Min(), i2da.J_Max(), 1,k_size, i2da.Halo_Size());
   for(int i=i_min_w_h;i<=i_max_w_h;i++) 
     for(int j=j_min_w_h;j<=j_max_w_h;j++)
       for(int k=k_min_w_h;k<=k_max_w_h;k++)
	 if(k==1) (*this)(i,j,1) = i2da(i,j); else (*this)(i,j,k) = (T)0;
  }

  ARRAY_3D(const ARRAY_3D<VECTOR_3D<T> >& va, const int component)
  {assert(component >=1 && component <=3); // component = 1 || 2 || 3
   Init_Array(va.I_Min(), va.I_Max(), va.J_Min(), va.J_Max(), 
	 	    va.K_Min(), va.K_Max(), va.Halo_Size());
   for(int n=0;n<total_size;n++) 
     switch(component){
       case 1: array[n] = va.Raw_Array(n).x; break;
       case 2: array[n] = va.Raw_Array(n).y; break;
       case 3: array[n] = va.Raw_Array(n).z; break;
     }
  } 

  ~ARRAY_3D() { delete [] array;}


  void Init_Array(const int imin,const int imax, const int jmin,const int jmax,
                  const int kmin,const int kmax, const int num_halo_cells);
  void Set_All_Elements_To(const T& element_value);
  void Set_Elements_To(const T* a, int size);
  void Set_NonHalo_Elements_To(const T& element_value);

  inline T& operator() (const int i, const int j, const int k);
  inline T  operator() (const int i, const int j, const int k) const;

  bool operator==(const ARRAY_3D<T>& a) const
  {if(Equal_Dimensions(*this,a)){
     for(int n=0;n<total_size;n++) if(array[n]!=a.Raw_Array(n)) return false;
     return true;}
   else return false;}

  ARRAY_3D<T>& operator=(const ARRAY_3D<T>& source)
  {assert(Equal_Dimensions(*this,source));
   for(int n=0;n<total_size;n++) array[n] = source.Raw_Array(n); return *this;}

  ARRAY_3D<T>& operator+=(const ARRAY_3D<T>& source)
  {assert(Equal_Dimensions(*this,source));
   for(int n=0;n<total_size;n++) array[n] += source.Raw_Array(n); return *this;}

  ARRAY_3D<T>& operator-=(const ARRAY_3D<T>& source)
  {assert(Equal_Dimensions(*this,source));
   for(int n=0;n<total_size;n++) array[n] -= source.Raw_Array(n); return *this;}

  ARRAY_3D<T>& operator*=(const ARRAY_3D<T>& source) //element-wise
  {assert(Equal_Dimensions(*this,source));
   for(int n=0;n<total_size;n++) array[n] *= source.Raw_Array(n); return *this;}

  ARRAY_3D<T>& operator*=(const T& scalar)
  {for(int n=0;n<total_size;n++) array[n] *= scalar; return *this;}

  ARRAY_3D<T> operator*(const T& scalar)
  {ARRAY_3D<T> a_copy(*this);
   for(int n=0;n<a_copy.total_size;n++) a_copy.Raw_Array(n) *= scalar; 
   return a_copy;}
  
  ARRAY_3D<T> operator*(const ARRAY_3D<T>& a)
  {assert(Equal_Dimensions(*this,a));
   ARRAY_3D<T> a_copy(*this);
   for(int n=0;n<total_size;n++) a_copy.Raw_Array(n) *= a.Raw_Array(n); 
   return a_copy;}

  ARRAY_3D<T> operator/(const ARRAY_3D<T>& a)
  {assert(Equal_Dimensions(*this,a));
   ARRAY_3D<T> a_copy(*this);
   for(int n=0;n<total_size;n++) {
     assert(a.Raw_Array(n)!=(T)0);
     a_copy.Raw_Array(n) /= a.Raw_Array(n);
   } 
   return a_copy;}

 /*
  friend ARRAY_3D<T> operator*(const T& scalar, ARRAY_3D<T>& a)
  {ARRAY_3D<T> a_copy(a);
   for(int n=0;n<a_copy.total_size;n++) a_copy.Raw_Array(n) *= scalar; 
   return a_copy;}
 */

  ARRAY_3D<T>& operator/=(const T scalar)
  {assert(scalar!=(T)0);
  for(int n=0;n<total_size;n++) array[n] /= scalar; return *this;}

  template<class T2>
  ARRAY_3D<T>& operator*=(const ARRAY_3D<T2>& av)
  {assert(Equal_Dimensions(*this,av));
   for(int n=0;n<total_size;n++) array[n] *= av.Raw_Array(n); return *this;}

  template<class T2>
  ARRAY_3D<T> operator*(const ARRAY_3D<T2>& av)
  {assert(Equal_Dimensions(*this,av));
   ARRAY_3D<T> a_copy(*this);
   for(int n=0;n<total_size;n++) a_copy.Raw_Array(n) *= av.Raw_Array(n); 
   return a_copy;}
  
  template<class T2>
  ARRAY_3D<T2> operator*(ARRAY_3D<T2>& av)
  {assert(Equal_Dimensions(*this,av));
   ARRAY_3D<T2> av_copy(av);   
   for(int n=0;n<total_size;n++) av_copy.Raw_Array(n)*=array[n]; 
   return av_copy;}

  template<class T2>
  friend ARRAY_3D<T> operator*(const T2 s, const ARRAY_3D<T>& a)
  {ARRAY_3D<T> a_copy(a);
   for(int n=0;n<a_copy.total_size;n++) a_copy.Raw_Array(n) *= s; 
   return a_copy;}

  ARRAY_3D<T> operator+(const ARRAY_3D<T>& a)
  {assert(Equal_Dimensions(*this,a));
   ARRAY_3D<T> a_copy(*this);
   for(int n=0;n<total_size;n++) a_copy.Raw_Array(n) += a.Raw_Array(n); 
   return a_copy;}

  template<class T2>
  static bool Equal_Dimensions(const ARRAY_3D<T>& a, const ARRAY_3D<T2>& b)
  {return a.I_Min()==b.I_Min() && a.I_Max()==b.I_Max() && a.J_Min()==b.J_Min()
       && a.J_Max()==b.J_Max() && a.K_Min()==b.K_Min() && a.K_Max()==b.K_Max()
       && a.Halo_Size()==b.Halo_Size();}

  int I_Min() const {return i_min;}
  int I_Max() const {return i_max;}
  int J_Min() const {return j_min;}
  int J_Max() const {return j_max;}
  int K_Min() const {return k_min;}
  int K_Max() const {return k_max;}
  int I_Min_With_Halo() const {return i_min_w_h;}
  int I_Max_With_Halo() const {return i_max_w_h;}
  int J_Min_With_Halo() const {return j_min_w_h;}
  int J_Max_With_Halo() const {return j_max_w_h;}
  int K_Min_With_Halo() const {return k_min_w_h;}
  int K_Max_With_Halo() const {return k_max_w_h;}

  int I_Size() const {return i_max - i_min + 1;}
  int J_Size() const {return j_max - j_min + 1;}
  int K_Size() const {return k_max - k_min + 1;}
  int Total_Size() const {return I_Size()*J_Size()*K_Size();}

  int I_Size_With_Halo() const {return i_size_w_h;}
  int J_Size_With_Halo() const {return j_size_w_h;}
  int K_Size_With_Halo() const {return k_size_w_h;}
  int Total_Size_With_Halo() const {return total_size;}
  int Halo_Size() const {return halo_size;}

  T& Raw_Array(const int n) {return array[n];}
  T Raw_Array(const int n) const {return array[n];}
  T* Raw_Array_Pointer() {return array;}
  void Delete_Array() {delete [] array;}

  bool Has_NAN_Values() 
  {for(int n=0;n<total_size;n++) if(isnan(array[n])) return true; return false;}

 private:
  int i_min, i_max, j_min, j_max, k_min, k_max; // w/o halo regions
  int i_min_w_h, i_max_w_h, j_min_w_h, j_max_w_h, k_min_w_h, k_max_w_h;
  int i_size_w_h, j_size_w_h, k_size_w_h, halo_size;
  int total_size;
  T*  array;
};

// global function
template <class T>
inline ostream& operator<< (ostream& output, const ARRAY_3D<T>& a)
{
  for (int i = a.I_Min_With_Halo(); i <= a.I_Max_With_Halo(); i++)
    for (int j = a.J_Min_With_Halo(); j <= a.J_Max_With_Halo(); j++)
      for (int k = a.K_Min_With_Halo(); k <= a.K_Max_With_Halo(); k++)
	output <<"["<<i<<","<<j<<","<<k<<"]"<<"="<<a(i,j,k)<<endl;
  return output;
}

template <class T>
inline T& ARRAY_3D<T>::operator () (const int i, const int j, const int k) 
{
  assert(i >= i_min_w_h && i <= i_max_w_h && j >= j_min_w_h && 
         j <= j_max_w_h && k >= k_min_w_h && k <= k_max_w_h);
  return array[i-i_min_w_h + (j-j_min_w_h)*i_size_w_h + 
		 (k-k_min_w_h)*i_size_w_h*j_size_w_h];
}

template <class T>
inline T ARRAY_3D<T>::operator () (const int i, const int j, const int k) const
{
  assert(i >= i_min_w_h && i <= i_max_w_h && j >= j_min_w_h && 
         j <= j_max_w_h && k >= k_min_w_h && k <= k_max_w_h);
  return array[i-i_min_w_h + (j-j_min_w_h)*i_size_w_h + 
		 (k-k_min_w_h)*i_size_w_h*j_size_w_h];
}
#endif
