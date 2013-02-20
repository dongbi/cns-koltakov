// vector_3d.h
// 3d vector class

#ifndef __VECTOR_3D__
#define __VECTOR_3D__

#include <iostream>
#include <assert.h>
using namespace std;

template<class T>
class VECTOR_3D
{
 public:
  VECTOR_3D() : x(T()), y(T()), z(T())
  {}

  VECTOR_3D(const T x_in, const T y_in, const T z_in)
            : x(x_in), y(y_in), z(z_in)
  {}

  VECTOR_3D(const T n_in)
            : x(n_in), y(n_in), z(n_in)
  {}

  VECTOR_3D(const VECTOR_3D<T>& v_in) : x(v_in.x), y(v_in.y), z(v_in.z)
  {}

  VECTOR_3D<T> operator-() const
  {return VECTOR_3D<T>(-x,-y,-z);}

  VECTOR_3D<T> operator*(const T scalar) const
  {return VECTOR_3D<T>(scalar*x, scalar*y, scalar*z);}

  VECTOR_3D<T> operator*(const VECTOR_3D<T>& v) const
  {return VECTOR_3D<T>(x*v.x, y*v.y, z*v.z);}

  VECTOR_3D<T>& operator*=(const T scalar)
  {x*=scalar; y*=scalar; z*=scalar; return *this;}

  VECTOR_3D<T>& operator*=(const VECTOR_3D<T>& v)
  {x*=v.x; y*=v.y; z*=v.z; return *this;}

  VECTOR_3D<T> operator/(const T scalar) const
  {return VECTOR_3D<T>(x/scalar, y/scalar, z/scalar);}

  VECTOR_3D<T>& operator/=(const T scalar)
  {assert(scalar!=(T)0); x/=scalar; y/=scalar; z/=scalar; return *this;}

  VECTOR_3D<T>& operator/=(const VECTOR_3D<T>& v)
  {x/=v.x; y/=v.y; z/=v.z; return *this;}

  friend VECTOR_3D<T> operator*(const T& scalar, const VECTOR_3D<T>& v)
  {return VECTOR_3D<T>(scalar*v.x, scalar*v.y, scalar*v.z);}

  VECTOR_3D<T> operator+(const VECTOR_3D<T>& v) const
  {return VECTOR_3D<T>(x+v.x,y+v.y,z+v.z);}

  VECTOR_3D<T> operator-(const VECTOR_3D<T>& v) const
  {return VECTOR_3D<T>(x-v.x,y-v.y,z-v.z);}

  VECTOR_3D<T>& operator+=(const VECTOR_3D<T>& v)
  {x+=v.x;y+=v.y; z+=v.z; return *this;}

  VECTOR_3D<T>& operator-=(const VECTOR_3D<T>& v)
  {x-=v.x; y-=v.y; z-=v.z; return *this;}

  bool operator== (const VECTOR_3D<T>& rhs) const
  {return x == rhs.x && y == rhs.y && z == rhs.z;}

  bool operator!= (const VECTOR_3D<T>& rhs) const
  {return x != rhs.x || y != rhs.y || z != rhs.z;}

  T Magnitude() const
  {return sqrt(x*x + y*y + z*z);}
  
  static T Dot_Product(const VECTOR_3D<T>& v1, const VECTOR_3D<T>& v2)
  {return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;}
 
  static VECTOR_3D<T> Cross_Product(
                                 const VECTOR_3D<T>& v1, const VECTOR_3D<T>& v2)
  {return VECTOR_3D<T>(v1.y*v2.z-v2.y*v1.z, 
                      -v1.x*v2.z+v2.x*v1.z, v1.x*v2.y-v2.x*v1.y);}

  T x,y,z;
};

template <class T>
inline bool isnan(const VECTOR_3D<T>& v)
{return isnan(v.x) || isnan(v.y) || isnan(v.z);}

template <class T>
inline ostream& operator<< (ostream& output, const VECTOR_3D<T>& v)
{output << "(" << v.x << "," << v.y << "," << v.z << ")"; return output;}

#endif
