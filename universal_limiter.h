// universal_limiter.h
// Universal Limiter classes to be used with TVD convection schemes
// in convection.h

#ifndef __UNIVERSAL_LIMITER__
#define __UNIVERSAL_LIMITER__

#include <cmath>
//*****************************************************************************
// Abstract class
//*****************************************************************************
template<class T=double>
class UNIVERSAL_LIMITER
{
 public:
  virtual T operator() (T ratio) const = 0;
};
//*****************************************************************************
template<class T=double>
class UPWIND_LIMITER : public UNIVERSAL_LIMITER<T>
{
 public:
  T operator() (T ratio) const {return (T)0;}
};
//*****************************************************************************
template<class T=double>
class MUSCL_LIMITER : public UNIVERSAL_LIMITER<T>
{
 public:
  T operator() (T ratio) const {
    return fmax( (T)0, fmin((T)2, fmin((T)2*ratio,(T).5*((T)1+ratio))) );
  }
};
//*****************************************************************************
#endif
