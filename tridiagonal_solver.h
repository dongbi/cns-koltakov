// tridiagonal_solver.h

#ifndef __TRIDIAGONAL_SOLVER__
#define __TRIDIAGONAL_SOLVER__

#include "mpi_driver.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "array_3d.h"
#include "array_2d.h"

#define PRINT_RESULTS

using namespace std;

template <class T=double> 
class TRIDIAGONAL_SOLVER
{
 public:
  TRIDIAGONAL_SOLVER(MPI_DRIVER<T>& md) 
                : mpi_comm(md.Get_Grid_Communicator()), my_rank(md.my_rank) {}
  ~TRIDIAGONAL_SOLVER() {}

  // wrapper for solvers with SCALAR RHS
  void Solve_Array_Of_Tridiagonal_Linear_Systems(ARRAY_2D<T>& a, ARRAY_2D<T>& b,ARRAY_2D<T>& c, ARRAY_2D<T>& f, bool periodic, int prev_neighbor, int next_neighbor)
   {
     if(periodic) Solve_Periodic_LS(a,b,c,f,prev_neighbor,next_neighbor);
     else Solve_Non_Periodic_LS(a,b,c,f,prev_neighbor,next_neighbor);
   }

  // wrapper for solvers with VECTOR_3D RHS
  void Solve_Array_Of_Tridiagonal_Linear_Systems(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c, ARRAY_2D<VECTOR_3D<T> >& f, bool periodic, int prev_neighbor, int next_neighbor)
  {
    if(periodic) Solve_Periodic_LS(a,b,c,f,prev_neighbor,next_neighbor);
    else Solve_Non_Periodic_LS(a,b,c,f,prev_neighbor,next_neighbor);
  }

 private:
  // non-periodic solvers: scalar and vector both use array[] version
  void Solve_Non_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c,
        ARRAY_2D<T>& f, int prev_neighbor, int next_neighbor);
  void Solve_Non_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c, 
        ARRAY_2D<VECTOR_3D<T> >& f, int prev_neighbor, int next_neighbor);
  void Solve_Non_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c,
        ARRAY_2D<T*>& f, int rhs_size, int prev_neighbor, int next_neighbor);

  // periodic solvers: scalar and vector both use array[] version
  void Solve_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c,
        ARRAY_2D<T>& f, int prev_neighbor, int next_neighbor);
  void Solve_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c, 
        ARRAY_2D<VECTOR_3D<T> >& f, int prev_neighbor, int next_neighbor);
  void Solve_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c,
        ARRAY_2D<T*>& f, int rhs_size, int prev_neighbor, int next_neighbor);

  MPI_Comm& mpi_comm;
  int my_rank;
};
//***************************************************************************
// LS solver: 1) Helper, 2) Non-periodic, 3) with ARRAY[] RHS (length 1,2,3,4)
//
// Note: 'f' is 2D array of arrays. Its length can be 1(scalar),3(vector) for 
// non-periodic case and 2(scalar+1),4(vector+1) for periodic case.
//****************************************************************************
template <class T>
void TRIDIAGONAL_SOLVER<T>::Solve_Non_Periodic_LS(ARRAY_2D<T>& a,ARRAY_2D<T>& b,
                                                ARRAY_2D<T>& c, ARRAY_2D<T*>& f,
                             int rhs_size, int prev_neighbor, int next_neighbor)
{
  assert(a.I_Min()==b.I_Min() && b.I_Min()==c.I_Min() && c.I_Min()==f.I_Min());
  assert(a.I_Max()==b.I_Max() && b.I_Max()==c.I_Max() && c.I_Max()==f.I_Max());
  assert(a.J_Min()==b.J_Min() && b.J_Min()==c.J_Min() && c.J_Min()==f.J_Min());
  assert(a.J_Max()==b.J_Max() && b.J_Max()==c.J_Max() && c.J_Max()==f.J_Max());

  int imin_w_h = a.I_Min_With_Halo(), imax_w_h = a.I_Max_With_Halo(), 
      jmin = a.J_Min_With_Halo()+1,   jmax = a.J_Max_With_Halo()-1,  
      isize_w_h = imax_w_h - imin_w_h + 1; 
  int size_of_border_msg        = isize_w_h * (rhs_size+1), 
      size_of_border_msg_short  = isize_w_h * rhs_size;
  double border_message[isize_w_h][rhs_size+1], 
         border_message_short[isize_w_h][rhs_size];
  MPI_Status status;

  if(prev_neighbor != MPI_PROC_NULL){
    MPI_Recv(border_message, size_of_border_msg, MPI_DOUBLE, 
	     prev_neighbor, 0, mpi_comm, &status);
    for(int i=imin_w_h;i<=imax_w_h;i++){
      for(int k=0;k<rhs_size;k++) f(i,jmin-1)[k] =border_message[i-imin_w_h][k];
      c(i,jmin-1) = border_message[i-imin_w_h][rhs_size];
    }
  }else // no previous neighbor
    for(int i=imin_w_h;i<=imax_w_h;i++){
      assert(b(i,jmin-1) != 0);
      for(int k=0;k<rhs_size;k++) f(i,jmin-1)[k] /= b(i,jmin-1);
      c(i,jmin-1) /= b(i,jmin-1);
    }

  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmin;j<=jmax;j++){
      b(i,j) -= a(i,j)*c(i,j-1);
      assert(b(i,j) != 0);
      b(i,j) = (T)1 / b(i,j);
      c(i,j) = c(i,j) * b(i,j);
    }

  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmin;j<=jmax;j++)
      for(int k=0;k<rhs_size;k++)
	f(i,j)[k] = (f(i,j)[k] - a(i,j)*f(i,j-1)[k]) * b(i,j);       

  if(next_neighbor != MPI_PROC_NULL){
    for(int i=imin_w_h;i<=imax_w_h;i++){
      for(int k=0;k<rhs_size;k++) border_message[i-imin_w_h][k] = f(i,jmax)[k];
      border_message[i-imin_w_h][rhs_size] = c(i,jmax);
    }
    MPI_Send(border_message, size_of_border_msg, MPI_DOUBLE, 
             next_neighbor, 0, mpi_comm);
  }else
    for(int i=imin_w_h;i<=imax_w_h;i++)
      for(int k=0;k<rhs_size;k++) {
	double denominator = b(i,jmax+1) - a(i,jmax+1)*c(i,jmax);
	assert(denominator != 0);
	f(i,jmax+1)[k] -=  a(i,jmax+1) * f(i,jmax)[k];
	f(i,jmax+1)[k] /= denominator;
      }
  
  if(next_neighbor != MPI_PROC_NULL){
    MPI_Recv(border_message_short, size_of_border_msg_short, MPI_DOUBLE, 
	     next_neighbor, 1, mpi_comm, &status);
     for(int i=imin_w_h;i<=imax_w_h;i++)
       for(int k=0;k<rhs_size;k++)
	 f(i,jmax+1)[k] = border_message_short[i-imin_w_h][k];
  }

  // backsolve
  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmax;j>=jmin;j--)
      for(int k=0;k<rhs_size;k++)
	f(i,j)[k] -= c(i,j)*f(i,j+1)[k];

  if(prev_neighbor != MPI_PROC_NULL){
    for(int i=imin_w_h;i<=imax_w_h;i++) 
      for(int k=0;k<rhs_size;k++)
	border_message_short[i-imin_w_h][k] = f(i,jmin)[k];
    MPI_Send(border_message_short, size_of_border_msg_short, MPI_DOUBLE, 
	     prev_neighbor, 1, mpi_comm);
  }else
    for(int i=imin_w_h;i<=imax_w_h;i++)
      for(int k=0;k<rhs_size;k++) 
	f(i,jmin-1)[k] -=  c(i,jmin-1) * f(i,jmin)[k];
}
//****************************************************************************
// Linear System solver: 1) Non-periodic, 2) with VECTOR_3D RHS 
//****************************************************************************
template <class T>
void TRIDIAGONAL_SOLVER<T>::Solve_Non_Periodic_LS(
   ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c, ARRAY_2D<VECTOR_3D<T> >& f,
   int prev_neighbor, int next_neighbor)
{
  // convert VECTOR_3D<T> -> T*: vector->array[0-2]
  ARRAY_2D<T*> f_array(f);
  for(int i=f_array.I_Min_With_Halo(); i<=f_array.I_Max_With_Halo(); i++)
    for(int j=f_array.J_Min_With_Halo(); j<=f_array.J_Max_With_Halo(); j++){
      f_array(i,j) = new T[3];
      f_array(i,j)[0] = f(i,j).x;
      f_array(i,j)[1] = f(i,j).y;
      f_array(i,j)[2] = f(i,j).z;
  }
  // call helper solver for array[] of size 3
  Solve_Non_Periodic_LS(a,b,c, f_array, 3, prev_neighbor,next_neighbor);
  // convert back T* -> VECTOR_3D<T>: array[0-2]->vector
  for(int i=f_array.I_Min_With_Halo(); i<=f_array.I_Max_With_Halo(); i++)
    for(int j=f_array.J_Min_With_Halo(); j<=f_array.J_Max_With_Halo(); j++){
      f(i,j).x = f_array(i,j)[0];
      f(i,j).y = f_array(i,j)[1];
      f(i,j).z = f_array(i,j)[2];
      delete [] f_array(i,j);
  }
  /*
  //DUPLICATION of code
  assert(a.I_Min()==b.I_Min() && b.I_Min()==c.I_Min() && c.I_Min()==f.I_Min());
  assert(a.I_Max()==b.I_Max() && b.I_Max()==c.I_Max() && c.I_Max()==f.I_Max());
  assert(a.J_Min()==b.J_Min() && b.J_Min()==c.J_Min() && c.J_Min()==f.J_Min());
  assert(a.J_Max()==b.J_Max() && b.J_Max()==c.J_Max() && c.J_Max()==f.J_Max());

  int imin_w_h = a.I_Min_With_Halo(),   imax_w_h = a.I_Max_With_Halo(), 
      jmin     = a.J_Min_With_Halo()+1, jmax     = a.J_Max_With_Halo()-1, 
      isize_w_h = imax_w_h-imin_w_h+1;

  int size_of_border_message  = isize_w_h * (3+1); //3d vector + 1
  double border_message[isize_w_h][3+1], border_message_short[isize_w_h][3];
  MPI_Status status;

  if(prev_neighbor != MPI_PROC_NULL){
    MPI_Recv(border_message, size_of_border_message, MPI_DOUBLE, 
	     prev_neighbor, 0, mpi_comm, &status);
    for(int i=imin_w_h;i<=imax_w_h;i++){
      f(i,jmin-1) = VECTOR_3D<T>(border_message[i-imin_w_h][0],
                 border_message[i-imin_w_h][1], border_message[i-imin_w_h][2]);
      c(i,jmin-1) = border_message[i-imin_w_h][3];
    }
  }else // no previous neighbor
    for(int i=imin_w_h;i<=imax_w_h;i++){
      assert(b(i,jmin-1) != 0);
      f(i,jmin-1) /= b(i,jmin-1);
      c(i,jmin-1) /= b(i,jmin-1);
    }

  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmin;j<=jmax;j++){
      b(i,j) -= a(i,j)*c(i,j-1);
      assert(b(i,j) != 0);
      b(i,j) = (T)1 / b(i,j);
      c(i,j) = c(i,j) * b(i,j);
    }

  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmin;j<=jmax;j++)
	f(i,j) = (f(i,j) - a(i,j)*f(i,j-1)) * b(i,j);       

  if(next_neighbor != MPI_PROC_NULL){
    size_of_border_message  = isize_w_h * (3+1);
    for(int i=imin_w_h;i<=imax_w_h;i++){
      border_message[i-imin_w_h][0] = f(i,jmax).x;
      border_message[i-imin_w_h][1] = f(i,jmax).y;
      border_message[i-imin_w_h][2] = f(i,jmax).z;
      border_message[i-imin_w_h][3] = c(i,jmax);
    }
    MPI_Send(border_message, size_of_border_message, MPI_DOUBLE, 
	     next_neighbor, 0, mpi_comm);
  }else
    for(int i=imin_w_h;i<=imax_w_h;i++){
	T denominator = b(i,jmax+1) - a(i,jmax+1)*c(i,jmax);
	assert(denominator != 0);
	f(i,jmax+1) -=  a(i,jmax+1) * f(i,jmax);
	f(i,jmax+1) /= denominator;
      }
  
  if(next_neighbor != MPI_PROC_NULL){
    size_of_border_message = isize_w_h*3;
    MPI_Recv(border_message_short, size_of_border_message, MPI_DOUBLE, 
	     next_neighbor, 1, mpi_comm, &status);
     for(int i=imin_w_h;i<=imax_w_h;i++)
       f(i,jmax+1) = VECTOR_3D<T>(border_message_short[i-imin_w_h][0],
                                  border_message_short[i-imin_w_h][1],
                                  border_message_short[i-imin_w_h][2]);
  }

  // backsolve
  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmax;j>=jmin;j--)
	f(i,j) -= c(i,j)*f(i,j+1);
  
  if(prev_neighbor != MPI_PROC_NULL){
    size_of_border_message = isize_w_h*3;
    for(int i=imin_w_h;i<=imax_w_h;i++) {
	border_message_short[i-imin_w_h][0] = f(i,jmin).x;
	border_message_short[i-imin_w_h][1] = f(i,jmin).y;	
	border_message_short[i-imin_w_h][2] = f(i,jmin).z;
    }
    MPI_Send(border_message_short, size_of_border_message, MPI_DOUBLE, 
	     prev_neighbor, 1, mpi_comm);
  }else
    for(int i=imin_w_h;i<=imax_w_h;i++)
	f(i,jmin-1) -=  c(i,jmin-1) * f(i,jmin);
  */
}
//****************************************************************************
// Linear System solver: 1) Non-periodic, 2) with SCALAR RHS 
//****************************************************************************
template <class T>
void TRIDIAGONAL_SOLVER<T>::Solve_Non_Periodic_LS(
   ARRAY_2D<T>& a, ARRAY_2D<T>& b, ARRAY_2D<T>& c, ARRAY_2D<T>& f,
   int prev_neighbor, int next_neighbor)
{
  // convert T->T*: scalar->array[0]
  ARRAY_2D<T*> f_array(f);
  for(int i=f_array.I_Min_With_Halo(); i<=f_array.I_Max_With_Halo(); i++)
    for(int j=f_array.J_Min_With_Halo(); j<=f_array.J_Max_With_Halo(); j++){
      f_array(i,j) = new T[1];
      f_array(i,j)[0] = f(i,j);
  }
  // call helper solver for array[] of size 1
  Solve_Non_Periodic_LS(a,b,c, f_array, 1, prev_neighbor,next_neighbor);
  // convert back T*->T: array[0]->scalar
  for(int i=f_array.I_Min_With_Halo(); i<=f_array.I_Max_With_Halo(); i++)
    for(int j=f_array.J_Min_With_Halo(); j<=f_array.J_Max_With_Halo(); j++){
      f(i,j) = f_array(i,j)[0];
      delete [] f_array(i,j);
  }
  /*
  //DUPLICATION of code
  assert(a.I_Min()==b.I_Min() && b.I_Min()==c.I_Min() && c.I_Min()==f.I_Min());
  assert(a.I_Max()==b.I_Max() && b.I_Max()==c.I_Max() && c.I_Max()==f.I_Max());
  assert(a.J_Min()==b.J_Min() && b.J_Min()==c.J_Min() && c.J_Min()==f.J_Min());
  assert(a.J_Max()==b.J_Max() && b.J_Max()==c.J_Max() && c.J_Max()==f.J_Max());
 
  int imin_w_h = a.I_Min_With_Halo(),   imax_w_h = a.I_Max_With_Halo(), 
      jmin     = a.J_Min_With_Halo()+1, jmax     = a.J_Max_With_Halo()-1, 
      isize_w_h = imax_w_h-imin_w_h+1;

  int size_of_border_message  = isize_w_h * (1+1); //scalar + 1
  double border_message[isize_w_h][1+1], border_message_short[isize_w_h];
  MPI_Status status;

  if(prev_neighbor != MPI_PROC_NULL){
    MPI_Recv(border_message, size_of_border_message, MPI_DOUBLE, 
	     prev_neighbor, 0, mpi_comm, &status);
    for(int i=imin_w_h;i<=imax_w_h;i++){
      f(i,jmin-1) = border_message[i-imin_w_h][0];
      c(i,jmin-1) = border_message[i-imin_w_h][1];
    }
  }else // no previous neighbor
    for(int i=imin_w_h;i<=imax_w_h;i++){
      assert(b(i,jmin-1) != 0);
      f(i,jmin-1) /= b(i,jmin-1);
      c(i,jmin-1) /= b(i,jmin-1);
    }

  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmin;j<=jmax;j++){
      b(i,j) -= a(i,j)*c(i,j-1);
      assert(b(i,j) != 0);
      b(i,j) = (T)1 / b(i,j);
      c(i,j) = c(i,j) * b(i,j);
    }

  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmin;j<=jmax;j++)
	f(i,j) = (f(i,j) - a(i,j)*f(i,j-1)) * b(i,j);       

  if(next_neighbor != MPI_PROC_NULL){
    size_of_border_message  = isize_w_h * (1+1);
    for(int i=imin_w_h;i<=imax_w_h;i++){
      border_message[i-imin_w_h][0] = f(i,jmax);
      border_message[i-imin_w_h][1] = c(i,jmax);
    }
    MPI_Send(border_message, size_of_border_message, MPI_DOUBLE, 
	     next_neighbor, 0, mpi_comm);
  }else
    for(int i=imin_w_h;i<=imax_w_h;i++){
	T denominator = b(i,jmax+1) - a(i,jmax+1)*c(i,jmax);
	assert(denominator != 0);
	f(i,jmax+1) -=  a(i,jmax+1) * f(i,jmax);
	f(i,jmax+1) /= denominator;
      }
  
  if(next_neighbor != MPI_PROC_NULL){
    size_of_border_message = isize_w_h; //short message
    MPI_Recv(border_message_short, size_of_border_message, MPI_DOUBLE, 
	     next_neighbor, 1, mpi_comm, &status);
     for(int i=imin_w_h;i<=imax_w_h;i++)
       f(i,jmax+1) = border_message_short[i-imin_w_h];
  }

  // backsolve
  for(int i=imin_w_h;i<=imax_w_h;i++)
    for(int j=jmax;j>=jmin;j--)
	f(i,j) -= c(i,j)*f(i,j+1);

  if(prev_neighbor != MPI_PROC_NULL){
    size_of_border_message = isize_w_h; //short message
    for(int i=imin_w_h;i<=imax_w_h;i++)
	border_message_short[i-imin_w_h] = f(i,jmin);
    MPI_Send(border_message_short, size_of_border_message, MPI_DOUBLE, 
	     prev_neighbor, 1, mpi_comm);
  }else
    for(int i=imin_w_h;i<=imax_w_h;i++)
	f(i,jmin-1) -=  c(i,jmin-1) * f(i,jmin);
  */
}
//****************************************************************************
// Linear System solver: 1) Periodic, 2) with 3D VECTOR RHS 
//****************************************************************************
template <class T>
void TRIDIAGONAL_SOLVER<T>::Solve_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, 
                                  ARRAY_2D<T>& c, ARRAY_2D<VECTOR_3D<T> >& f,
                                          int prev_neighbor, int next_neighbor)
{
  int prev_neighbor_mod = prev_neighbor, next_neighbor_mod = next_neighbor,
      imin_w_h  = a.I_Min_With_Halo(), imax_w_h = a.I_Max_With_Halo(), 
      jmin_w_h  = a.J_Min_With_Halo(), jmax_w_h = a.J_Max_With_Halo();
  // for every element of 2D array: 1) converting vector -> array[]
  //                                2) adding an extra element (RHS)
  ARRAY_2D<T*> f_mod(f);
  for(int i=imin_w_h; i<=imax_w_h; i++)
    for(int j=jmin_w_h; j<=jmax_w_h; j++){
      f_mod(i,j) = new T[4];
      f_mod(i,j)[0] = f(i,j).x;
      f_mod(i,j)[1] = f(i,j).y;
      f_mod(i,j)[2] = f(i,j).z;
      f_mod(i,j)[3] = (T)0;
  }

  // call helper solver for array[] of size 3+1
  Solve_Periodic_LS(a,b,c, f_mod, 3, prev_neighbor_mod, next_neighbor_mod);

  // for each element of f_mod array: 1) copy array[0] -> scalar
  //                                  2) delete array[]
  for(int i=imin_w_h; i<=imax_w_h; i++) 
    for(int j=jmin_w_h; j<=jmax_w_h; j++){
      f(i,j).x = f_mod(i,j)[0];
      f(i,j).y = f_mod(i,j)[1];
      f(i,j).z = f_mod(i,j)[2];
      delete [] f_mod(i,j);
  }
}
//****************************************************************************
// Linear System solver: 1) Periodic, 2) with SCALAR RHS 
//****************************************************************************
template <class T>
void TRIDIAGONAL_SOLVER<T>::Solve_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, 
                                              ARRAY_2D<T>& c, ARRAY_2D<T>& f, 
                                          int prev_neighbor, int next_neighbor)
{
  int prev_neighbor_mod = prev_neighbor, next_neighbor_mod = next_neighbor,
      imin_w_h  = a.I_Min_With_Halo(), imax_w_h = a.I_Max_With_Halo(), 
      jmin_w_h  = a.J_Min_With_Halo(), jmax_w_h = a.J_Max_With_Halo();
  // for every element of 2D array: 1) copy scalar -> array[0]
  //                                2) adding an extra element (RHS)
  ARRAY_2D<T*> f_mod(f);
  for(int i=imin_w_h; i<=imax_w_h; i++)
    for(int j=jmin_w_h; j<=jmax_w_h; j++){
      f_mod(i,j) = new T[2];
      f_mod(i,j)[0] = f(i,j);
      f_mod(i,j)[1] = (T)0;
  }

  // call helper solver for array[] of size 1+1
  Solve_Periodic_LS(a,b,c, f_mod, 1, prev_neighbor_mod, next_neighbor_mod);

  // for each element of f_mod array: 1) copy array[0] -> scalar
  //                                  2) delete array[]
  for(int i=imin_w_h; i<=imax_w_h; i++) 
    for(int j=jmin_w_h; j<=jmax_w_h; j++){
      f(i,j) = f_mod(i,j)[0];
      delete [] f_mod(i,j);
  }
}
//****************************************************************************
// LS solver: 1) Helper, 2) Periodic, 3) with ARRAY[] RHS (length 2 or 4)
//****************************************************************************
template <class T>
void TRIDIAGONAL_SOLVER<T>::Solve_Periodic_LS(ARRAY_2D<T>& a, ARRAY_2D<T>& b, 
                                              ARRAY_2D<T>& c, ARRAY_2D<T*>& f, 
			     int rhs_size, int prev_neighbor, int next_neighbor)
{
  int prev_neighbor_mod = prev_neighbor, next_neighbor_mod = next_neighbor,
      imin_w_h  = a.I_Min_With_Halo(), imax_w_h = a.I_Max_With_Halo(), 
      jmin_w_h  = a.J_Min_With_Halo(), jmax_w_h = a.J_Max_With_Halo(), 
      isize_w_h = imax_w_h - imin_w_h + 1;
  double send_msg[isize_w_h][rhs_size+1], send_msg_short[isize_w_h][rhs_size],
         recv_msg[isize_w_h][rhs_size+1], recv_msg_short[isize_w_h][rhs_size];
  MPI_Status status;

  if(prev_neighbor >= my_rank){
    prev_neighbor_mod = MPI_PROC_NULL;
    for(int i=imin_w_h; i<=imax_w_h; i++){
      for(int k=0; k<rhs_size+1; k++) f(i,jmin_w_h)[k] = (T)0; //k=1..L+1
      a(i,jmin_w_h) = (T)0;
      b(i,jmin_w_h) = (T)1;
      c(i,jmin_w_h) = (T)0;
      f(i,jmin_w_h+1)[rhs_size] = a(i,jmin_w_h+1); //extra element (L+1)
      b(i,jmin_w_h+1) -= a(i,jmin_w_h+1);
      a(i,jmin_w_h+1) = (T)0;
    }
  }

  if(next_neighbor <= my_rank){
    next_neighbor_mod = MPI_PROC_NULL;
    for(int i=imin_w_h; i<=imax_w_h; i++){
      for(int k=0; k<rhs_size+1; k++) f(i,jmax_w_h)[k] = (T)0; //k=1..L+1
      a(i,jmax_w_h) = (T)0;
      b(i,jmax_w_h) = (T)1;
      c(i,jmax_w_h) = (T)0;
      f(i,jmax_w_h-1)[rhs_size] = c(i,jmax_w_h-1); //extra element (L+1)
      b(i,jmax_w_h-1) -= c(i,jmax_w_h-1);
      c(i,jmax_w_h-1) = (T)0;
    }
  }

  Solve_Non_Periodic_LS(a,b,c,f,rhs_size+1,prev_neighbor_mod,next_neighbor_mod);

  if(prev_neighbor == my_rank){
    for(int i=imin_w_h; i<=imax_w_h; i++)
      for(int k=0; k<rhs_size+1; k++)  
	recv_msg[i-imin_w_h][k] = f(i,jmax_w_h-1)[k] + f(i,jmin_w_h+1)[k];
    for(int i=imin_w_h; i<=imax_w_h; i++){ 
      T denominator = (T)1 + recv_msg[i-imin_w_h][rhs_size];
      //cout<<"Denom1:"<<denominator<<"fm="<<f(i,jmax_w_h-1)[rhs_size]<<"f1="<<f(i,jmin_w_h+1)[rhs_size]<<endl;
      assert(denominator);
      for(int k=0; k<rhs_size; k++)
	send_msg_short[i-imin_w_h][k] = recv_msg[i-imin_w_h][k] / denominator;
     }
  }

  if(prev_neighbor > my_rank){
    MPI_Recv(recv_msg, isize_w_h*(rhs_size+1), MPI_DOUBLE, 
	     prev_neighbor, 0, mpi_comm, &status);    
    for(int i=imin_w_h; i<=imax_w_h; i++){
      for(int k=0; k<rhs_size+1; k++)
	recv_msg[i-imin_w_h][k] += f(i,jmin_w_h+1)[k];
      for(int k=0; k<rhs_size; k++){
        T denominator = (T)1 + recv_msg[i-imin_w_h][rhs_size];     
        //cout<<"Denom2:"<<denominator<<endl;
        assert(denominator);
	send_msg_short[i-imin_w_h][k] = recv_msg[i-imin_w_h][k] / denominator;
      }
    }
  }

  if(next_neighbor < my_rank){
    for(int i=imin_w_h; i<=imax_w_h; i++)
      for(int k=0; k<rhs_size+1; k++)    
	send_msg[i-imin_w_h][k] = f(i,jmax_w_h-1)[k];
    MPI_Send(send_msg, isize_w_h*(rhs_size+1), MPI_DOUBLE, 
	     next_neighbor, 0, mpi_comm);   
  }

  if(prev_neighbor < my_rank){
    MPI_Recv(recv_msg_short, isize_w_h*rhs_size, MPI_DOUBLE, 
	     prev_neighbor, 0, mpi_comm, &status);     
    for(int i=imin_w_h; i<=imax_w_h; i++)
      for(int k=0; k<rhs_size; k++)
	send_msg_short[i-imin_w_h][k] = recv_msg_short[i-imin_w_h][k];
  }

  if(next_neighbor > my_rank)
    MPI_Send(send_msg_short, isize_w_h*rhs_size, MPI_DOUBLE, 
	     next_neighbor, 0, mpi_comm);      
  
  for(int i=imin_w_h; i<=imax_w_h; i++)
    for(int j=jmin_w_h; j<=jmax_w_h; j++)    
      for(int k=0; k<rhs_size; k++)    
	f(i,j)[k] -= send_msg_short[i-imin_w_h][k] * f(i,j)[rhs_size];
}
//****************************************************************************
#endif
