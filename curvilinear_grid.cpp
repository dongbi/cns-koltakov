//*****************************************************************************
// File:  curvilinear_grid.cpp
// Class: CURVILINEAR_GRID<T>
// Note:  Curvilinear grid class
//*****************************************************************************
#include "curvilinear_grid.h"
#include <cmath>
#include <assert.h>
//*****************************************************************************
  template<class T>
  CURVILINEAR_GRID<T>::CURVILINEAR_GRID(PARAMETERS<T> *params, MPI_DRIVER<T> *mpi)
: mpi_driver(mpi), parameters(params)
{
  Init_Empty_Grid(params->x_min, params->x_max,params->y_min, 
      params->y_max, params->z_min, params->z_max, params->i_min,params->i_max,
      params->j_min, params->j_max, params->k_min, params->k_max,
      params->halo_size, params->mg_sub_levels);

  //local grid: [0..1; 0..1; 0..1] or read from file
  if(!params->read_grid_from_file){
    Init_Local_Grid_Node_Positions_With_Uniform_Cube();
    Custom_Adjust_Local_Grid_Node_Positions();
  }else 
    mpi_driver->Read_Global_Array_From_Disk(params->grid_filename, *grid);
  //Init_Local_Grid_Node_Positions_From_File(); //serial version

  Finish_Initialization_Based_On_Node_Positions();
}
//*****************************************************************************
// Finish static grid initialization: invoked by constructor
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Finish_Initialization_Based_On_Node_Positions()
{
  cout<<"Proc#"<<mpi_driver->my_rank<<": Calculating Metrics..."<<endl;
  Calculate_Metrics();

  cout<<"Proc#"<<mpi_driver->my_rank<<": Creating Subgrids..."<<endl;
  Create_Update_Subgrids();
} 
//*****************************************************************************
// Class variables initialization: invoked by constructor
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Init_Empty_Grid(T xmin,T xmax, T ymin,T ymax, T zmin,
    T zmax, int imin, int imax, int jmin, int jmax, int kmin, int kmax,
    int nhalo, int subgrid_levs)
{
  x_min=xmin; x_max=xmax; y_min=ymin; y_max=ymax; z_min=zmin; z_max=zmax;
  i_size = imax-imin+1; j_size = jmax-jmin+1; k_size = kmax-kmin+1; 

  x_length = xmax-xmin; y_length = ymax-ymin; z_length = zmax-zmin;
  subgrid_levels = subgrid_levs;

  halo_size = nhalo; i_size_w_h = i_size+2*nhalo; 
  j_size_w_h = j_size+2*nhalo; k_size_w_h = k_size+2*nhalo;

  i_min=imin; i_max=imax; j_min=jmin; j_max=jmax; k_min=kmin; k_max=kmax;
  i_min_w_h = i_min-halo_size; i_max_w_h = i_max+halo_size;
  j_min_w_h = j_min-halo_size; j_max_w_h = j_max+halo_size;
  k_min_w_h = k_min-halo_size; k_max_w_h = k_max+halo_size;

  grid = 
    new ARRAY_3D<VECTOR_3D<T> >(i_min,i_max,j_min,j_max,k_min,k_max,halo_size);
  //allocating metric structures
  XI_x = new ARRAY_3D<T>(*grid); XI_y = new ARRAY_3D<T>(*grid);
  XI_z = new ARRAY_3D<T>(*grid);
  ET_x = new ARRAY_3D<T>(*grid); ET_y = new ARRAY_3D<T>(*grid);
  ET_z = new ARRAY_3D<T>(*grid);
  ZT_x = new ARRAY_3D<T>(*grid); ZT_y = new ARRAY_3D<T>(*grid);
  ZT_z = new ARRAY_3D<T>(*grid);
  inverse_Jacobian = new ARRAY_3D<T>(*grid);
  G11 = new ARRAY_3D<T>(*grid); G12 = new ARRAY_3D<T>(*grid);
  G13 = new ARRAY_3D<T>(*grid);
  G21 = new ARRAY_3D<T>(*grid); G22 = new ARRAY_3D<T>(*grid);
  G23 = new ARRAY_3D<T>(*grid);
  G31 = new ARRAY_3D<T>(*grid); G32 = new ARRAY_3D<T>(*grid);
  G33 = new ARRAY_3D<T>(*grid); GCC = new ARRAY_3D<T>(*grid);
}
//*****************************************************************************
// Destructor
//*****************************************************************************
  template<class T> 
CURVILINEAR_GRID<T>::~CURVILINEAR_GRID()
{
  delete grid; 
  delete XI_x; delete XI_y; delete XI_z;
  delete ET_x; delete ET_y; delete ET_z;
  delete ZT_x; delete ZT_y; delete ZT_z;
  delete inverse_Jacobian;
  delete G11; delete G12; delete G13;  delete G21; delete G22; delete G23;
  delete G31; delete G32; delete G33;  delete GCC;
  if(subgrid_levels){
    delete num_x_sub; delete num_y_sub; delete num_z_sub;
    for (int level = 1; level <= subgrid_levels; level++) {
      delete (*inv_Jac_sub)(level);
      delete (*G11_sub)(level); delete (*G12_sub)(level);
      delete (*G13_sub)(level); delete (*G21_sub)(level); 
      delete (*G22_sub)(level); delete (*G23_sub)(level);
      delete (*G31_sub)(level); delete (*G32_sub)(level);
      delete (*G33_sub)(level); delete (*GCC_sub)(level);
    }
    delete inv_Jac_sub; delete GCC_sub; 
    delete G11_sub; delete G12_sub; delete G13_sub;
    delete G21_sub; delete G22_sub; delete G23_sub;
    delete G31_sub; delete G32_sub; delete G33_sub;
  }
}
//*****************************************************************************
// Inits node positions for local grid given owning Proc's Cartesian coords
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Init_Local_Grid_Node_Positions_With_Uniform_Cube()
{  
  // local domain min
  int i_local_min = mpi_driver->my_coords_in_grid[0] * i_size,
      j_local_min = mpi_driver->my_coords_in_grid[1] * j_size, 
      k_local_min = mpi_driver->my_coords_in_grid[2] * k_size;

  //Make global nodes cover 3d cube: [0..1; 0..1; 0..1]
  for(int i = i_min_w_h; i <= i_max_w_h; i++)
    for(int j = j_min_w_h; j <= j_max_w_h; j++)
      for(int k = k_min_w_h; k <= k_max_w_h; k++){
        (*grid)(i,j,k).x = 
          (T(i_local_min+i) - .5) / (T)parameters->num_total_nodes_x;
        (*grid)(i,j,k).y = 
          (T(j_local_min+j) - .5) / (T)parameters->num_total_nodes_y;
        (*grid)(i,j,k).z = 
          (T(k_local_min+k) - .5) / (T)parameters->num_total_nodes_z;
      }
}
//*****************************************************************************
// Stretches the horizontal (x) coordinate to resolve breaking: assumes x:[0;1]
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Stretch_In_Horizontal_To_Resolve_Breaking()
{
  int i_local_min = mpi_driver->my_coords_in_grid[0] * i_size;
  T r = 1. - (parameters->x_stretching_ratio - 1.), 
  dx_max = (1. - r) / (1. - pow(r,parameters->num_total_nodes_x-1));
  for(int i = i_min_w_h; i <= i_max_w_h; i++)
    for(int j = j_min_w_h; j <= j_max_w_h; j++)
      for(int k = k_min_w_h; k <= k_max_w_h; k++)
        if(i_local_min+i-1 >= i_min_w_h)
          (*grid)(i,j,k).x = dx_max * (1. - pow(r,i_local_min+i-1)) / (1. - r);
        else
          (*grid)(i,j,k).x = -dx_max * halo_size; //halo equal cells below 0
}
//*****************************************************************************
// Stretches the vertical (z) coordinate to resolve the bottom: assumes y:[0;1]
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Stretch_In_Vertical_To_Resolve_Bottom()
{
  int k_local_min = mpi_driver->my_coords_in_grid[2] * k_size;
  T r = parameters->z_stretching_ratio, 
  dz_min = (1. - r) / (1. - pow(r,parameters->num_total_nodes_z-1));
  for(int i = i_min_w_h; i <= i_max_w_h; i++)
    for(int j = j_min_w_h; j <= j_max_w_h; j++)
      for(int k = k_min_w_h; k <= k_max_w_h; k++)
        if(k_local_min+k-1 >= k_min_w_h)
          (*grid)(i,j,k).z = dz_min * (1. - pow(r,k_local_min+k-1)) / (1. - r);
        else
          (*grid)(i,j,k).z = -dz_min * halo_size; //halo equal cells below 0
}
//*****************************************************************************
// Reads node locations from file: 1 Proc version
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Init_Local_Grid_Node_Positions_From_File()
{
  if(!mpi_driver->my_rank){
    ifstream infile("grid.dat", ios::in|ios::binary);
    if(infile.is_open()){
      int size_i, size_j, size_k;
      infile.read((char*)(&size_i), sizeof(int)); assert(size_i==i_size_w_h);
      infile.read((char*)(&size_j), sizeof(int)); assert(size_j==j_size_w_h);
      infile.read((char*)(&size_k), sizeof(int)); assert(size_k==k_size_w_h);

      for(int i = i_min_w_h; i <= i_max_w_h; i++)
        for(int j = j_min_w_h; j <= j_max_w_h; j++)
          for(int k = k_min_w_h; k <= k_max_w_h; k++){
            infile.read((char*)(&(*grid)(i,j,k).x),sizeof(double));
            infile.read((char*)(&(*grid)(i,j,k).y),sizeof(double));
            infile.read((char*)(&(*grid)(i,j,k).z),sizeof(double));
          }
      infile.close();
    }else cout<<"ERROR: Cannot open the grid file"<<endl;
  }
}
//*****************************************************************************
// Adjust the cube to fit the domain along with any requested streching
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Custom_Adjust_Local_Grid_Node_Positions()
{
  //Shift the domain in case the skew angle is not 90 degrees
  if(parameters->domain_skew_angle < (T)90) 
    Shift_Grid_Into_Parallelogram(parameters->domain_skew_angle);

  //Increase resolution at the boundaries
  if(parameters->stretch_in_x || parameters->stretch_in_y
      || parameters->stretch_in_z) 
    Push_Nodes_Toward_Boundaries();

  if(parameters->x_stretching_ratio) Stretch_In_Horizontal_To_Resolve_Breaking();
  if(parameters->z_stretching_ratio) Stretch_In_Vertical_To_Resolve_Bottom();

  //Set depth(:,:) based on node positions (modified from uniform grid)
  parameters->Set_Depth_Based_On_Node_Locations(*this->grid);

  //Make global nodes cover 3d cube: [xmin..xmax; ymin..ymax; zmin..zmax]
  //or for variable depth: [xmin..xmax; depth(:,:)..0; zmin..zmax]
  Scale_And_Shift_Grid_Nodes_To_Fit_Physical_Domain();

  //Increase resolution around an interface in Y (push nodes towards the middle)
  if(parameters->resolve_interface_in_z){
    ARRAY_1D<T> interface_in_Z(i_size_w_h, i_min_w_h,(T)1);
    Adjust_Nodes_To_Resolve_Interface_In_Z(interface_in_Z);
  }
}
//*****************************************************************************
// Calculation of Inverse Jacobian and Mesh Skewness Tensor
// given node coordinates.
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Calculate_Metrics()
{
  T X_xi,Y_xi,Z_xi, X_et,Y_et,Z_et, X_zt,Y_zt,Z_zt, inv_Jac, 
    XI_s_X,ET_s_X,ZT_s_X, XI_s_Y,ET_s_Y,ZT_s_Y, XI_s_Z,ET_s_Z,ZT_s_Z;

  // I-face calculations
  for (int i=i_min_w_h; i<i_max_w_h; i++)
    for (int j=j_min_w_h + 1; j<j_max_w_h; j++)
      for (int k=k_min_w_h + 1; k<k_max_w_h; k++){
        X_xi = (*grid)(i+1,j,k).x - (*grid)(i,j,k).x;
        Y_xi = (*grid)(i+1,j,k).y - (*grid)(i,j,k).y;
        Z_xi = (*grid)(i+1,j,k).z - (*grid)(i,j,k).z;

        X_et = (T).25 * ( (*grid)(i,  j+1,k).x - (*grid)(i,  j-1,k).x + 
            (*grid)(i+1,j+1,k).x - (*grid)(i+1,j-1,k).x );
        Y_et = (T).25 * ( (*grid)(i,  j+1,k).y - (*grid)(i,  j-1,k).y + 
            (*grid)(i+1,j+1,k).y - (*grid)(i+1,j-1,k).y );
        Z_et = (T).25 * ( (*grid)(i,  j+1,k).z - (*grid)(i,  j-1,k).z + 
            (*grid)(i+1,j+1,k).z - (*grid)(i+1,j-1,k).z );

        X_zt = (T).25 * ( (*grid)(i,  j,k+1).x - (*grid)(i,  j,k-1).x + 
            (*grid)(i+1,j,k+1).x - (*grid)(i+1,j,k-1).x );
        Y_zt = (T).25 * ( (*grid)(i,  j,k+1).y - (*grid)(i,  j,k-1).y + 
            (*grid)(i+1,j,k+1).y - (*grid)(i+1,j,k-1).y );
        Z_zt = (T).25 * ( (*grid)(i,  j,k+1).z - (*grid)(i,  j,k-1).z + 
            (*grid)(i+1,j,k+1).z - (*grid)(i+1,j,k-1).z );

        inv_Jac =  (X_xi*Y_et*Z_zt + X_et*Y_zt*Z_xi + X_zt*Y_xi*Z_et - 
            X_zt*Y_et*Z_xi - X_et*Y_xi*Z_zt - X_xi*Y_zt*Z_et );
        /*if(inv_Jac == 0)
          cout<<"I = "<<i<<",J = "<<j<<",K = "<<k
          <<",X_xi="<<X_xi<<",Y_et="<<Y_et<<",Z_zt="<<Z_zt<<endl;*/
        assert(inv_Jac);
        inv_Jac = (T)1 / inv_Jac;

        XI_s_X = ( Y_et * Z_zt - Y_zt * Z_et );
        ET_s_X = ( Y_zt * Z_xi - Y_xi * Z_zt );
        ZT_s_X = ( Y_xi * Z_et - Y_et * Z_xi );
        XI_s_Y = ( Z_et * X_zt - Z_zt * X_et );
        ET_s_Y = ( Z_zt * X_xi - Z_xi * X_zt );
        ZT_s_Y = ( Z_xi * X_et - Z_et * X_xi );
        XI_s_Z = ( X_et * Y_zt - X_zt * Y_et );
        ET_s_Z = ( X_zt * Y_xi - X_xi * Y_zt );
        ZT_s_Z = ( X_xi * Y_et - X_et * Y_xi );

        (*G11)(i,j,k) = inv_Jac*(XI_s_X*XI_s_X + XI_s_Y*XI_s_Y + XI_s_Z*XI_s_Z);
        (*G12)(i,j,k) = inv_Jac*(XI_s_X*ET_s_X + XI_s_Y*ET_s_Y + XI_s_Z*ET_s_Z);
        (*G13)(i,j,k) = inv_Jac*(XI_s_X*ZT_s_X + XI_s_Y*ZT_s_Y + XI_s_Z*ZT_s_Z);

        (*XI_x)(i,j,k) = XI_s_X;
        (*XI_y)(i,j,k) = XI_s_Y;
        (*XI_z)(i,j,k) = XI_s_Z;
      }

  // J-face calculations
  for (int i=i_min_w_h + 1; i<i_max_w_h; i++)
    for (int j=j_min_w_h; j<j_max_w_h; j++)
      for (int k=k_min_w_h + 1; k<k_max_w_h; k++){
        X_et = (*grid)(i,j+1,k).x - (*grid)(i,j,k).x;
        Y_et = (*grid)(i,j+1,k).y - (*grid)(i,j,k).y;
        Z_et = (*grid)(i,j+1,k).z - (*grid)(i,j,k).z;

        X_zt = (T).25 * ( (*grid)(i,j,  k+1).x - (*grid)(i,j,  k-1).x + 
            (*grid)(i,j+1,k+1).x - (*grid)(i,j+1,k-1).x );
        Y_zt = (T).25 * ( (*grid)(i,j,  k+1).y - (*grid)(i,j,  k-1).y + 
            (*grid)(i,j+1,k+1).y - (*grid)(i,j+1,k-1).y );
        Z_zt = (T).25 * ( (*grid)(i,j,  k+1).z - (*grid)(i,j,  k-1).z + 
            (*grid)(i,j+1,k+1).z - (*grid)(i,j+1,k-1).z );

        X_xi = (T).25 * ( (*grid)(i+1,j,  k).x - (*grid)(i-1,j,  k).x + 
            (*grid)(i+1,j+1,k).x - (*grid)(i-1,j+1,k).x );
        Y_xi = (T).25 * ( (*grid)(i+1,j,  k).y - (*grid)(i-1,j,  k).y + 
            (*grid)(i+1,j+1,k).y - (*grid)(i-1,j+1,k).y );
        Z_xi = (T).25 * ( (*grid)(i+1,j,  k).z - (*grid)(i-1,j,  k).z + 
            (*grid)(i+1,j+1,k).z - (*grid)(i-1,j+1,k).z );

        inv_Jac = (X_xi*Y_et*Z_zt + X_et*Y_zt*Z_xi + X_zt*Y_xi*Z_et - 
            X_zt*Y_et*Z_xi - X_et*Y_xi*Z_zt - X_xi*Y_zt*Z_et);
        assert(inv_Jac);
        inv_Jac = (T)1 / inv_Jac;

        XI_s_X = ( Y_et * Z_zt - Y_zt * Z_et );
        ET_s_X = ( Y_zt * Z_xi - Y_xi * Z_zt );
        ZT_s_X = ( Y_xi * Z_et - Y_et * Z_xi );
        XI_s_Y = ( Z_et * X_zt - Z_zt * X_et );
        ET_s_Y = ( Z_zt * X_xi - Z_xi * X_zt );
        ZT_s_Y = ( Z_xi * X_et - Z_et * X_xi );
        XI_s_Z = ( X_et * Y_zt - X_zt * Y_et );
        ET_s_Z = ( X_zt * Y_xi - X_xi * Y_zt );
        ZT_s_Z = ( X_xi * Y_et - X_et * Y_xi );

        (*G21)(i,j,k) = inv_Jac*(ET_s_X*XI_s_X + ET_s_Y*XI_s_Y + ET_s_Z*XI_s_Z);
        (*G22)(i,j,k) = inv_Jac*(ET_s_X*ET_s_X + ET_s_Y*ET_s_Y + ET_s_Z*ET_s_Z);
        (*G23)(i,j,k) = inv_Jac*(ET_s_X*ZT_s_X + ET_s_Y*ZT_s_Y + ET_s_Z*ZT_s_Z);

        (*ET_x)(i,j,k) = ET_s_X;
        (*ET_y)(i,j,k) = ET_s_Y;
        (*ET_z)(i,j,k) = ET_s_Z;
      }

  // K-face calculations
  for (int i=i_min_w_h + 1; i<i_max_w_h; i++)
    for (int j=j_min_w_h + 1; j<j_max_w_h; j++)
      for (int k=k_min_w_h; k<k_max_w_h; k++){
        X_zt = (*grid)(i,j,k+1).x - (*grid)(i,j,k).x;
        Y_zt = (*grid)(i,j,k+1).y - (*grid)(i,j,k).y;
        Z_zt = (*grid)(i,j,k+1).z - (*grid)(i,j,k).z;

        X_xi = (T).25 * ( (*grid)(i+1,j,k  ).x - (*grid)(i-1,j,k  ).x + 
            (*grid)(i+1,j,k+1).x - (*grid)(i-1,j,k+1).x );
        Y_xi = (T).25 * ( (*grid)(i+1,j,k  ).y - (*grid)(i-1,j,k  ).y + 
            (*grid)(i+1,j,k+1).y - (*grid)(i-1,j,k+1).y );
        Z_xi = (T).25 * ( (*grid)(i+1,j,k  ).z - (*grid)(i-1,j,k  ).z + 
            (*grid)(i+1,j,k+1).z - (*grid)(i-1,j,k+1).z );

        X_et = (T).25 * ( (*grid)(i,j+1,k  ).x - (*grid)(i,j-1,k  ).x + 
            (*grid)(i,j+1,k+1).x - (*grid)(i,j-1,k+1).x );
        Y_et = (T).25 * ( (*grid)(i,j+1,k  ).y - (*grid)(i,j-1,k  ).y + 
            (*grid)(i,j+1,k+1).y - (*grid)(i,j-1,k+1).y );
        Z_et = (T).25 * ( (*grid)(i,j+1,k  ).z - (*grid)(i,j-1,k  ).z + 
            (*grid)(i,j+1,k+1).z - (*grid)(i,j-1,k+1).z );

        inv_Jac = (X_xi*Y_et*Z_zt + X_et*Y_zt*Z_xi + X_zt*Y_xi*Z_et - 
            X_zt*Y_et*Z_xi - X_et*Y_xi*Z_zt - X_xi*Y_zt*Z_et);
        assert(inv_Jac);
        inv_Jac = (T)1 / inv_Jac;

        XI_s_X = ( Y_et * Z_zt - Y_zt * Z_et );
        ET_s_X = ( Y_zt * Z_xi - Y_xi * Z_zt );
        ZT_s_X = ( Y_xi * Z_et - Y_et * Z_xi );
        XI_s_Y = ( Z_et * X_zt - Z_zt * X_et );
        ET_s_Y = ( Z_zt * X_xi - Z_xi * X_zt );
        ZT_s_Y = ( Z_xi * X_et - Z_et * X_xi );
        XI_s_Z = ( X_et * Y_zt - X_zt * Y_et );
        ET_s_Z = ( X_zt * Y_xi - X_xi * Y_zt );
        ZT_s_Z = ( X_xi * Y_et - X_et * Y_xi );

        (*G31)(i,j,k) = inv_Jac*(ZT_s_X*XI_s_X + ZT_s_Y*XI_s_Y + ZT_s_Z*XI_s_Z);
        (*G32)(i,j,k) = inv_Jac*(ZT_s_X*ET_s_X + ZT_s_Y*ET_s_Y + ZT_s_Z*ET_s_Z);
        (*G33)(i,j,k) = inv_Jac*(ZT_s_X*ZT_s_X + ZT_s_Y*ZT_s_Y + ZT_s_Z*ZT_s_Z);

        (*ZT_x)(i,j,k) = ZT_s_X;
        (*ZT_y)(i,j,k) = ZT_s_Y;
        (*ZT_z)(i,j,k) = ZT_s_Z;
      }

  // Center calculations

  for (int i=i_min_w_h + 1; i<i_max_w_h; i++)
    for (int j=j_min_w_h + 1; j<j_max_w_h; j++)
      for (int k=k_min_w_h + 1; k<k_max_w_h; k++){
        X_xi = (T).5 * ( (*grid)(i+1,j,k).x - (*grid)(i-1,j,k).x );
        Y_xi = (T).5 * ( (*grid)(i+1,j,k).y - (*grid)(i-1,j,k).y );
        Z_xi = (T).5 * ( (*grid)(i+1,j,k).z - (*grid)(i-1,j,k).z );

        X_et = (T).5 * ( (*grid)(i,j+1,k).x - (*grid)(i,j-1,k).x );
        Y_et = (T).5 * ( (*grid)(i,j+1,k).y - (*grid)(i,j-1,k).y );
        Z_et = (T).5 * ( (*grid)(i,j+1,k).z - (*grid)(i,j-1,k).z );

        X_zt = (T).5 * ( (*grid)(i,j,k+1).x - (*grid)(i,j,k-1).x );
        Y_zt = (T).5 * ( (*grid)(i,j,k+1).y - (*grid)(i,j,k-1).y );
        Z_zt = (T).5 * ( (*grid)(i,j,k+1).z - (*grid)(i,j,k-1).z );

        inv_Jac = (X_xi*Y_et*Z_zt + X_et*Y_zt*Z_xi + X_zt*Y_xi*Z_et -
            X_zt*Y_et*Z_xi - X_et*Y_xi*Z_zt - X_xi*Y_zt*Z_et);	
        assert(inv_Jac);

        // save change in Jacobians (inv_Jac is not inverted just yet)
        //if(parameters->moving_grid && (*inverse_Jacobian)(i,j,k))
        //  (*Jacobian_diff)(i,j,k) = (inv_Jac-(T)1/(*inverse_Jacobian)(i,j,k))
        //                           / parameters->delta_time;

        (*inverse_Jacobian)(i,j,k) = (T)1 / inv_Jac;
      }

  for (int i=i_min_w_h; i<=i_max_w_h; i++)
    for (int j=j_min_w_h; j<=j_max_w_h; j++)
      for (int k=k_min_w_h; k<=k_max_w_h; k++){
        (*G12)(i,j,k) = (T).25 * (*G12)(i,j,k);
        (*G13)(i,j,k) = (T).25 * (*G13)(i,j,k);
        (*G23)(i,j,k) = (T).25 * (*G23)(i,j,k);
        (*G21)(i,j,k) = (T).25 * (*G21)(i,j,k);
        (*G31)(i,j,k) = (T).25 * (*G31)(i,j,k);
        (*G32)(i,j,k) = (T).25 * (*G32)(i,j,k);
        if( (i > i_min_w_h) && (j > j_min_w_h) && (k > k_min_w_h) )
          (*GCC)(i,j,k) = (*G11)(i,j,k) + (*G11)(i-1,j,k) +
            (*G22)(i,j,k) + (*G22)(i,j-1,k) +
            (*G33)(i,j,k) + (*G33)(i,j,k-1);
      }
}
//*****************************************************************************
// Creating (or updating) data structures for subgrids 
// with # of levels = subgrid_levels. Used for multigrid calculations.
// If 'create_data_structures' = false, structures are reused.
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Create_Update_Subgrids(bool create_data_structures)
{
  // quit if there are no multigrid sublevels
  if(!subgrid_levels){
    num_x_sub = num_y_sub = num_z_sub = NULL;
    inv_Jac_sub = G11_sub = G12_sub=G13_sub = G21_sub = G22_sub = G23_sub = 
      G31_sub = G32_sub = G33_sub = GCC_sub = NULL;
    return;
  }
  // create data structures
  if(create_data_structures){
    num_x_sub = new ARRAY_1D<int>(subgrid_levels);
    num_y_sub = new ARRAY_1D<int>(subgrid_levels);
    num_z_sub = new ARRAY_1D<int>(subgrid_levels);
    inv_Jac_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G11_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G12_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G13_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G21_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G22_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G23_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G31_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G32_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    G33_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
    GCC_sub = new ARRAY_1D<ARRAY_3D<T>* >(subgrid_levels);
  }

  // traverse these structres down to the lowest sub-level
  for (int level = 1; level <= subgrid_levels; level++) {
    // create data structres at sublevels
    if(create_data_structures){
      int nx,ny,nz, nx_sub,ny_sub,nz_sub;
      if (level == 1) {
        nx = i_size;
        ny = j_size;
        nz = k_size; 
      } else {
        nx = (*num_x_sub)(level-1);
        ny = (*num_y_sub)(level-1);
        nz = (*num_z_sub)(level-1);
      }

      nx_sub = nx / 2;
      nz_sub = nz / 2;
      if(!parameters->two_d)
        ny_sub = ny / 2;
      else
        ny_sub = ny;

      (*num_x_sub)(level) = nx_sub;
      (*num_y_sub)(level) = ny_sub;
      (*num_z_sub)(level) = nz_sub;

      (*G11_sub)(level) = new ARRAY_3D<T>(1,nx_sub,1,ny_sub,1,nz_sub,halo_size);
      (*G12_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G13_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G21_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G22_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G23_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G31_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G32_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*G33_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*GCC_sub)(level) = new ARRAY_3D<T>(*(*G11_sub)(level)); 
      (*inv_Jac_sub)(level)=new ARRAY_3D<T>(*(*G11_sub)(level)); 
    }

    // initialize data structres at sublevels
    if (level == 1)
      Subsample_Metric_Quantities ((*num_x_sub)(level), (*num_y_sub)(level),
          (*num_z_sub)(level),
          *G11,*G12,*G13,*G21,*G22,*G23,*G31,*G32,*G33,*GCC,*inverse_Jacobian,
          *(*G11_sub)(level),*(*G12_sub)(level),*(*G13_sub)(level),
          *(*G21_sub)(level),*(*G22_sub)(level),*(*G23_sub)(level),
          *(*G31_sub)(level),*(*G32_sub)(level),*(*G33_sub)(level),
          *(*GCC_sub)(level),*(*inv_Jac_sub)(level) );
    else
      Subsample_Metric_Quantities ((*num_x_sub)(level), (*num_y_sub)(level),
          (*num_z_sub)(level),
          *(*G11_sub)(level-1),*(*G12_sub)(level-1),*(*G13_sub)(level-1),
          *(*G21_sub)(level-1),*(*G22_sub)(level-1),*(*G23_sub)(level-1),
          *(*G31_sub)(level-1),*(*G32_sub)(level-1),*(*G33_sub)(level-1),
          *(*GCC_sub)(level-1),*(*inv_Jac_sub)(level-1),
          *(*G11_sub)(level),*(*G12_sub)(level),*(*G13_sub)(level),
          *(*G21_sub)(level),*(*G22_sub)(level),*(*G23_sub)(level),
          *(*G31_sub)(level),*(*G32_sub)(level),*(*G33_sub)(level),
          *(*GCC_sub)(level),*(*inv_Jac_sub)(level) );
  }
}
//*****************************************************************************
// Subsampling of the Inverse Jacobian and Mesh Skewness Tensor
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Subsample_Metric_Quantities(
    const int nx_new, const int ny_new, const int nz_new,
    const ARRAY_3D<T>& G11, const ARRAY_3D<T>& G12, const ARRAY_3D<T>& G13,
    const ARRAY_3D<T>& G21, const ARRAY_3D<T>& G22, const ARRAY_3D<T>& G23,
    const ARRAY_3D<T>& G31, const ARRAY_3D<T>& G32, const ARRAY_3D<T>& G33, 
    const ARRAY_3D<T>& GCC, const ARRAY_3D<T>& inverse_jacobian, 
    ARRAY_3D<T>& G11_new, ARRAY_3D<T>& G12_new, ARRAY_3D<T>& G13_new,
    ARRAY_3D<T>& G21_new, ARRAY_3D<T>& G22_new, ARRAY_3D<T>& G23_new,
    ARRAY_3D<T>& G31_new, ARRAY_3D<T>& G32_new, ARRAY_3D<T>& G33_new, 
    ARRAY_3D<T>& GCC_new, ARRAY_3D<T>& inverse_jacobian_new)

{
  int i_new, j_new, k_new;

  if(!parameters->two_d){
    for (int i=1; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          i_new = 2*i - 1; 
          j_new = 2*j - 1; 
          k_new = 2*k - 1;
          inverse_jacobian_new(i,j,k) = 
            1.0 / inverse_jacobian(i_new,  j_new,  k_new  )
          + 1.0 / inverse_jacobian(i_new+1,j_new,  k_new  )
          + 1.0 / inverse_jacobian(i_new,  j_new+1,k_new  )
          + 1.0 / inverse_jacobian(i_new+1,j_new+1,k_new  )
          + 1.0 / inverse_jacobian(i_new,  j_new,  k_new+1)
          + 1.0 / inverse_jacobian(i_new+1,j_new,  k_new+1)
          + 1.0 / inverse_jacobian(i_new,  j_new+1,k_new+1)
          + 1.0 / inverse_jacobian(i_new+1,j_new+1,k_new+1);
          inverse_jacobian_new(i,j,k) = 1.0 / inverse_jacobian_new(i,j,k);
        }
    for (int i=0; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          i_new = 2*i; 
          j_new = 2*j - 1; 
          k_new = 2*k - 1;
          G11_new(i,j,k) = 0.5 * ( G11(i_new,j_new,  k_new  )
                                 + G11(i_new,j_new+1,k_new  )
                                 + G11(i_new,j_new  ,k_new+1)
                                 + G11(i_new,j_new+1,k_new+1) );
          G12_new(i,j,k) = 0.5 * ( G12(i_new,j_new,  k_new  )
                                 + G12(i_new,j_new+1,k_new  )
                                 + G12(i_new,j_new  ,k_new+1)
                                 + G12(i_new,j_new+1,k_new+1) );
          G13_new(i,j,k) = 0.5 * ( G13(i_new,j_new,  k_new  )
                                 + G13(i_new,j_new+1,k_new  )
                                 + G13(i_new,j_new  ,k_new+1)
                                 + G13(i_new,j_new+1,k_new+1) );
        }
    for (int i=1; i<=nx_new; i++)
      for (int j=0; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          i_new = 2*i - 1; 
          j_new = 2*j; 
          k_new = 2*k - 1;
          G21_new(i,j,k) = 0.5 * ( G21(i_new,  j_new,k_new  )
                                 + G21(i_new+1,j_new,k_new  )
                                 + G21(i_new,  j_new,k_new+1)
                                 + G21(i_new+1,j_new,k_new+1) );
          G22_new(i,j,k) = 0.5 * ( G22(i_new,  j_new,k_new  )
                                 + G22(i_new+1,j_new,k_new  )
                                 + G22(i_new,  j_new,k_new+1)
                                 + G22(i_new+1,j_new,k_new+1) );
          G23_new(i,j,k) = 0.5 * ( G23(i_new,  j_new,k_new  )
                                 + G23(i_new+1,j_new,k_new  )
                                 + G23(i_new,  j_new,k_new+1)
                                 + G23(i_new+1,j_new,k_new+1) );
        }
    for (int i=1; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=0; k<=nz_new; k++) {
          i_new = 2*i - 1; 
          j_new = 2*j - 1; 
          k_new = 2*k;
          G31_new(i,j,k) = 0.5 * ( G31(i_new,  j_new,  k_new)
                                 + G31(i_new+1,j_new,  k_new)
                                 + G31(i_new,  j_new+1,k_new)
                                 + G31(i_new+1,j_new+1,k_new) );
          G32_new(i,j,k) = 0.5 * ( G32(i_new,  j_new,  k_new)
                                 + G32(i_new+1,j_new,  k_new)
                                 + G32(i_new,  j_new+1,k_new)
                                 + G32(i_new+1,j_new+1,k_new) );
          G33_new(i,j,k) = 0.5 * ( G33(i_new,  j_new,  k_new)
                                 + G33(i_new+1,j_new,  k_new)
                                 + G33(i_new,  j_new+1,k_new)
                                 + G33(i_new+1,j_new+1,k_new) );
        }
    for (int i=1; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          GCC_new(i,j,k) = G11_new(i,j,k) +  G11_new(i-1,j,  k  )
            + G22_new(i,j,k) +  G22_new(i,  j-1,k  )
            + G33_new(i,j,k) +  G33_new(i,  j,  k-1);
        }
  } else {
    for (int i=1; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          i_new = 2*i - 1; 
          j_new = j; 
          k_new = 2*k - 1;
          inverse_jacobian_new(i,j,k) = 
            1.0 / inverse_jacobian(i_new,  j_new,  k_new  )
          + 1.0 / inverse_jacobian(i_new+1,j_new,  k_new  )
          + 1.0 / inverse_jacobian(i_new,  j_new,  k_new+1)
          + 1.0 / inverse_jacobian(i_new+1,j_new,  k_new+1);
          inverse_jacobian_new(i,j,k) = 1.0 / inverse_jacobian_new(i,j,k);
        }
    for (int i=0; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          i_new = 2*i; 
          j_new = j; 
          k_new = 2*k - 1;
          G11_new(i,j,k) = 0.5 * ( G11(i_new,j_new,  k_new  )
                                 + G11(i_new,j_new  ,k_new+1) );
          G12_new(i,j,k) = 0.5 * ( G12(i_new,j_new,  k_new  )
                                 + G12(i_new,j_new  ,k_new+1) );
          G13_new(i,j,k) = 0.5 * ( G13(i_new,j_new,  k_new  )
                                 + G13(i_new,j_new  ,k_new+1) );
        }
    for (int i=1; i<=nx_new; i++)
      for (int j=0; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          i_new = 2*i - 1; 
          j_new = j; 
          k_new = 2*k - 1;
          G21_new(i,j,k) = 0.5 * ( G21(i_new,  j_new,k_new  )
                                 + G21(i_new+1,j_new,k_new  )
                                 + G21(i_new,  j_new,k_new+1)
                                 + G21(i_new+1,j_new,k_new+1) );
          G22_new(i,j,k) = 0.5 * ( G22(i_new,  j_new,k_new  )
                                 + G22(i_new+1,j_new,k_new  )
                                 + G22(i_new,  j_new,k_new+1)
                                 + G22(i_new+1,j_new,k_new+1) );
          G23_new(i,j,k) = 0.5 * ( G23(i_new,  j_new,k_new  )
                                 + G23(i_new+1,j_new,k_new  )
                                 + G23(i_new,  j_new,k_new+1)
                                 + G23(i_new+1,j_new,k_new+1) );
        }
    for (int i=1; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=0; k<=nz_new; k++) {
          i_new = 2*i - 1; 
          j_new = j; 
          k_new = 2*k;
          G31_new(i,j,k) = 0.5 * ( G31(i_new,  j_new,  k_new)
                                 + G31(i_new+1,j_new,  k_new) );
          G32_new(i,j,k) = 0.5 * ( G32(i_new,  j_new,  k_new)
                                 + G32(i_new+1,j_new,  k_new) );
          G33_new(i,j,k) = 0.5 * ( G33(i_new,  j_new,  k_new)
                                 + G33(i_new+1,j_new,  k_new) );
        }
    for (int i=1; i<=nx_new; i++)
      for (int j=1; j<=ny_new; j++)
        for (int k=1; k<=nz_new; k++) {
          GCC_new(i,j,k) = G11_new(i,j,k) +  G11_new(i-1,j,  k  )
            + G22_new(i,j,k) +  G22_new(i,  j-1,k  )
            + G33_new(i,j,k) +  G33_new(i,  j,  k-1);
        }
    }
}
//*****************************************************************************
// Calculate total volume of domain
//*****************************************************************************
  template<class T>
T CURVILINEAR_GRID<T>::Calculate_Domain_Volume()
{
  T total_volume = (T)0;
  for (int i=i_min; i<=i_max; i++)
    for (int j=j_min; j<=j_max; j++)
      for (int k=k_min; k<=k_max; k++)
        total_volume += (T)1 / (*inverse_Jacobian)(i,j,k);
  //sum all subdomains
  mpi_driver->Replace_With_Sum_On_All_Procs(total_volume);
  return total_volume;
}
//*****************************************************************************
// Clusters grid points around the interface in the vertical dim (depth)
// Location of the interface is determined by interface_Y(i),
// it is (= 1) for the undisturbed interface.
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Adjust_Nodes_To_Resolve_Interface_In_Z(
    const ARRAY_1D<T>& interface_Z)
{
  assert(interface_Z.Min_Index() == i_min_w_h);
  assert(interface_Z.Max_Index() == i_max_w_h);
  T ap = (T)4;
  T dp = 1.;//abs(z_min);
  for (int i=i_min_w_h; i<=i_max_w_h; i++)
    for (int j=j_min_w_h; j<=j_max_w_h; j++)
      for (int k=k_min_w_h; k<=k_max_w_h; k++){
        T z = (*grid)(i,j,k).z;
        (*grid)(i,j,k).z = interface_Z(i)
          * (exp(ap*z)-exp(-ap*dp)) / ((T)1-exp(-ap*dp));
        (*grid)(i,j,k).z -= ((T)2-interface_Z(i))
          * (exp(ap*(-1.-z))-exp(-ap*dp)) / ((T)1-exp(-ap*dp)) 
          + interface_Z(i);
        (*grid)(i,j,k).z *= .5*dp;
      }
}  
//*****************************************************************************
// Does not change the number of nodes but 
// alters their existing position in a uniform Cartesian cube [0,1;0,1;0,1].
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Push_Nodes_Toward_Boundaries()
{
  T am = (T)3.2,
    bm = (T).48,
    dm = (T)0,
    cm = (T)1/( dm  + log(cosh(am*(1 - bm)) / cosh(am*bm))     / am
        - log(cosh(am*bm)       / cosh(am*(1-bm))) / am ); 

  for (int i=i_min_w_h; i<=i_max_w_h; i++)
    for (int j=j_min_w_h; j<=j_max_w_h; j++)
      for (int k=k_min_w_h; k<=k_max_w_h; k++){  
        if(parameters->stretch_in_x)
          (*grid)(i,j,k).x = cm*( dm * (*grid)(i,j,k).x + 
              log( cosh(am*((*grid)(i,j,k).x - bm))   / cosh(am*bm) )    / am 
              - log( cosh(am*((*grid)(i,j,k).x - 1+bm)) / cosh(am*(1-bm))) / am );
        if(parameters->stretch_in_y)
          (*grid)(i,j,k).y = cm*( dm * (*grid)(i,j,k).y + 
              log( cosh(am*((*grid)(i,j,k).y - bm))   / cosh(am*bm) )    / am 
              - log( cosh(am*((*grid)(i,j,k).y - 1+bm)) / cosh(am*(1-bm))) / am );
        if(parameters->stretch_in_z)
          (*grid)(i,j,k).z = cm*( dm * (*grid)(i,j,k).z + 
              log( cosh(am*((*grid)(i,j,k).z - bm))   / cosh(am*bm) )    / am 
              - log( cosh(am*((*grid)(i,j,k).z - 1+bm)) / cosh(am*(1-bm))) / am );
      }
}
//*****************************************************************************
// Grid factory function
//*****************************************************************************
  template<class T>
ARRAY_3D<VECTOR_3D<T> >* CURVILINEAR_GRID<T>::Create_Uniform_Grid()
{
  T dx = x_length / (parameters->num_total_nodes_x-1),
    dy = y_length / (parameters->num_total_nodes_y-1),
    dz = z_length / (parameters->num_total_nodes_z-1);

  ARRAY_3D<VECTOR_3D<T> >* grid = 
    new ARRAY_3D<VECTOR_3D<T> >(i_min,i_max,j_min,j_max,k_min,k_max,halo_size);

  for (int i=i_min_w_h; i<=i_max_w_h; i++)
    for (int j=j_min_w_h; j<=j_max_w_h; j++)
      for (int k=k_min_w_h; k<=k_max_w_h; k++)
        (*grid)(i,j,k) = VECTOR_3D<T>(x_min+i*dx, y_min+j*dy, z_min+k*dz);
  return grid;
}
//*****************************************************************************
// Adjusts grid nodes according to the size of physical domain
// In case of variable depth, the y-coordinate is within [depth(:,:), 0]
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Scale_And_Shift_Grid_Nodes_To_Fit_Physical_Domain()
{
  for (int i=i_min_w_h; i<=i_max_w_h; i++)
    for (int j=j_min_w_h; j<=j_max_w_h; j++)
      for (int k=k_min_w_h; k<=k_max_w_h; k++){
        (*grid)(i,j,k).x *= x_length; 
        (*grid)(i,j,k).x += x_min; 
        if(parameters->depth){
          (*grid)(i,j,k).z *= (*parameters->depth)(i,j); 
          (*grid)(i,j,k).z -= (*parameters->depth)(i,j); 
        }else{
          (*grid)(i,j,k).z *= z_length; 
          (*grid)(i,j,k).z += z_min;
        }
        (*grid)(i,j,k).y *= y_length; 
        (*grid)(i,j,k).y += y_min;
      }
} 
//*****************************************************************************
// Creates parallelogram by shifting nodes in X with no variability in Y or Z
//*****************************************************************************
  template<class T>
void CURVILINEAR_GRID<T>::Shift_Grid_Into_Parallelogram(const T angle_in_deg)
{
  T dz = z_length / (parameters->num_total_nodes_z-1),
    shift_value = dz/tan((T)3.1415*angle_in_deg/(T)180);

  for(int i=i_min_w_h; i<=i_max_w_h; i++)
    for(int j=j_min_w_h; j<=j_max_w_h; j++)
      for(int k=k_min_w_h; k<=k_max_w_h; k++)
        (*grid)(i,j,k).x += shift_value*(k-k_min_w_h);

  // update x_max with the delta shift of the top layer
  parameters->x_max += shift_value*(k_max_w_h-k_min_w_h);
}
//*****************************************************************************
template class CURVILINEAR_GRID<double>;
