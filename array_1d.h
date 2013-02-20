//****************************************************************************
// File:    array_1d.h
// Class:   ARRAY_1D<T>
// Comment: 1D array class with index range: (min_n .. max_n=(min_n+size-1))
//****************************************************************************
#ifndef __ARRAY_1D__
#define __ARRAY_1D__
#include <iostream>
#include <assert.h>
using namespace std;

template <class T=double> 
class ARRAY_1D
{
 public:
 ARRAY_1D(const int size, const int min = 1, const T value = T());
  ~ARRAY_1D() { delete [] array;}

  void Set_All_Elements_To(const T element_value);
  int Min_Index() const {return min_n;}
  int Max_Index() const {return max_n;}
  int Size() const {return size;}
  T& operator() (const int n);
  T  operator() (const int n) const;

 private:
  int min_n, max_n, size;
  T*  array;
};
//****************************************************************************
template <class T>
ARRAY_1D<T>::ARRAY_1D(const int size, const int min, const T value) 
                     : size(size), min_n(min), max_n(min_n+size-1)
{
  array = new T[size]; 
  for(int n=0; n<size; n++) array[n] = value;
}
//***************************************************************************
template <class T>
void ARRAY_1D<T>::Set_All_Elements_To(const T element_value)
{
  for(int n=0; n<size; n++) array[n] = element_value;
}
//****************************************************************************
template <class T>
T& ARRAY_1D<T>::operator () (const int n) 
{
  assert(n >= min_n && n <= max_n);
  return array[n-min_n];
}
//****************************************************************************
template <class T>
T ARRAY_1D<T>::operator () (const int n) const 
{
  assert(n >= min_n && n <= max_n);
  return array[n-min_n];
}
//****************************************************************************
// global function
template <class T>
inline ostream& operator<< (ostream& output, const ARRAY_1D<T>& theArray)
{
  for (int n = theArray.Min_Index(); n <= theArray.Max_Index(); n++)
	output <<"["<<n<<"]"<<"="<<theArray(n)<<endl;
  return output;
}
//****************************************************************************
#endif
