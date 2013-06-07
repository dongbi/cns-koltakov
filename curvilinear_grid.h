//*****************************************************************************
// File:  curvilinear_grid.h
// Class: CURVILINEAR_GRID<T>
// Note:  Curvilinear grid in 3D
//*****************************************************************************
#ifndef __CURVILINEAR_GRID__
#define __CURVILINEAR_GRID__

#include "mpi_driver.h"
#include "array_3d.h"
#include "array_2d.h"
#include "array_1d.h"
#include "vector_3d.h"
#include "parameters.h"

//template<class T> class MPI_DRIVER;

template<class T=double>
class CURVILINEAR_GRID
{
 public:
  //ctor1: sets node positions to (0,0,0)
  CURVILINEAR_GRID(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax, int imin, 
    int imax,int jmin,int jmax,int kmin,int kmax, int nhalo, int subgrid_levels)
  {
    Init_Empty_Grid(xmin, xmax, ymin, ymax, zmin, zmax,
                 imin, imax, jmin, jmax, kmin, kmax, nhalo, subgrid_levels);
    parameters=NULL; mpi_driver=NULL;
  }
  //ctor2: creates distributed grid among procs in Cartesian topology
  CURVILINEAR_GRID(PARAMETERS<T> *parameters, MPI_DRIVER<T> *mpi_driver);
  ~CURVILINEAR_GRID();

  void Init_Empty_Grid(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax, int imin,
     int imax,int jmin,int jmax,int kmin,int kmax, int nhalo, int subgrid_levs);
  void Finish_Initialization_Based_On_Node_Positions();

  void Init_Local_Grid_Node_Positions_From_File();
  void Init_Local_Grid_Node_Positions_With_Uniform_Cube();
  void Custom_Adjust_Local_Grid_Node_Positions();

  void Calculate_Metrics();
  void Create_Update_Subgrids(bool create_data_structures=true);

  T Calculate_Domain_Volume();

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
  int I_Size() const {return i_size;} //nodes range:{i/j/k}_min..{i/j/k}_max
  int J_Size() const {return j_size;}
  int K_Size() const {return k_size;}
  int Halo_Size() const {return halo_size;}
  VECTOR_3D<T>& operator() (int i, int j, int k) {return (*grid)(i,j,k);}
  VECTOR_3D<T>  operator() (int i, int j, int k) const {return (*grid)(i,j,k);}

  ARRAY_3D<VECTOR_3D<T> >* grid;
  ARRAY_3D<T> *XI_x,*XI_y,*XI_z, *ET_x,*ET_y,*ET_z, *ZT_x,*ZT_y,*ZT_z,
              *XI_x_c,*XI_y_c,*XI_z_c, *ET_x_c,*ET_y_c,*ET_z_c, *ZT_x_c,*ZT_y_c,*ZT_z_c,
              *inverse_Jacobian,
              *G11,*G12,*G13, *G21,*G22,*G23, *G31,*G32,*G33, *GCC;

  // data structures for subgrids
  ARRAY_1D<int> *num_x_sub, *num_y_sub, *num_z_sub;
  ARRAY_1D<ARRAY_3D<T>* > *inv_Jac_sub, *G11_sub,*G12_sub,*G13_sub, 
                          *G21_sub,*G22_sub,*G23_sub, 
                          *G31_sub,*G32_sub,*G33_sub, *GCC_sub;
 private:
  void Stretch_In_Horizontal_To_Resolve_Breaking();
  void Stretch_In_Vertical_To_Resolve_Bottom();
  void Push_Nodes_Toward_Boundaries();
  void Scale_And_Shift_Grid_Nodes_To_Fit_Physical_Domain();
  void Adjust_Nodes_To_Resolve_Interface_In_Z(const ARRAY_1D<T>& interface_Z);
  void Shift_Grid_Into_Parallelogram(const T slope);

  void Subsample_Metric_Quantities(
        const int nx_new, const int ny_new, const int nz_new,
	const ARRAY_3D<T>& G11, const ARRAY_3D<T>& G12, const ARRAY_3D<T>& G13,
	const ARRAY_3D<T>& G21, const ARRAY_3D<T>& G22, const ARRAY_3D<T>& G23,
	const ARRAY_3D<T>& G31, const ARRAY_3D<T>& G32, const ARRAY_3D<T>& G33,
	const ARRAY_3D<T>& GCC, const ARRAY_3D<T>& jacobian, 
	      ARRAY_3D<T>& G11_new, ARRAY_3D<T>& G12_new, ARRAY_3D<T>& G13_new,
	      ARRAY_3D<T>& G21_new, ARRAY_3D<T>& G22_new, ARRAY_3D<T>& G23_new,
	      ARRAY_3D<T>& G31_new, ARRAY_3D<T>& G32_new, ARRAY_3D<T>& G33_new,
	      ARRAY_3D<T>& GCC_new, ARRAY_3D<T>& jacobian_new);

  ARRAY_3D<VECTOR_3D<T> >* Create_Uniform_Grid();

  T x_min, x_max, y_min, y_max, z_min, z_max, // global domain boundaries
    x_length, y_length, z_length; 
  // grid dimensions
  int i_size, j_size, k_size, i_size_w_h, j_size_w_h, k_size_w_h, halo_size,    
      i_min, i_max, j_min, j_max, k_min, k_max,
      i_min_w_h, i_max_w_h, j_min_w_h, j_max_w_h, k_min_w_h, k_max_w_h,//w/ halo
      subgrid_levels;

 protected:
  PARAMETERS<T> *parameters;
  MPI_DRIVER<T> *mpi_driver; // for distributed grids
};
//*****************************************************************************
#endif
