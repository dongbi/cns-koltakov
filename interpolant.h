//*****************************************************************************
// File:  interpolant.h
// Class: INTERPOLANT<T>
// Note:  global grid interpolation between grids
//*****************************************************************************
#ifndef __INTERPOLANT__
#define __INTERPOLANT__

#include "mpi_driver.h"

template<class T=double>
class INTERPOLANT
{
 public:
  INTERPOLANT(MPI_DRIVER<T> *mpi, 
              ARRAY_3D<VECTOR_3D<T> >& local_grid, ARRAY_3D<T>& local_f);
  ~INTERPOLANT();

  void Interpolate_On_New_Grid(ARRAY_3D<VECTOR_3D<T> >& new_grid, 
                               ARRAY_3D<T>& new_f);
 private:
  void Find_Confining_Cell_On_Global_Grid(int k, VECTOR_3D<T>& location, 
			        VECTOR_3D<int>& node1, VECTOR_3D<int>& node2);
  T Interpolate_Function_In_Quad(VECTOR_3D<T>& P,
			         VECTOR_3D<int>& node1, VECTOR_3D<int>& node2);
  T Interpolate_Function_In_Triangle(VECTOR_3D<T>& p,
	VECTOR_3D<T>& v1, VECTOR_3D<T>& v2, VECTOR_3D<T>& v3, T f1, T f2, T f3);
  MPI_DRIVER<T> *mpi_driver;
  ARRAY_3D<T> *global_F;
  ARRAY_3D<VECTOR_3D<T> > *global_grid;
};
//*****************************************************************************
// Constructor
//*****************************************************************************
template<class T>
INTERPOLANT<T>::INTERPOLANT(MPI_DRIVER<T> *mpi,
                     ARRAY_3D<VECTOR_3D<T> >& local_grid, ARRAY_3D<T>& local_F) 
  : mpi_driver(mpi)
{
  int local_x_size = local_grid.I_Max() - local_grid.I_Min() + 1, 
      local_y_size = local_grid.J_Max() - local_grid.J_Min() + 1,
      local_z_size = local_grid.K_Max() - local_grid.K_Min() + 1,
      halo_size = local_grid.Halo_Size(),
      local_x_size_w_h = local_x_size+2*halo_size,
      local_y_size_w_h = local_y_size+2*halo_size,
      local_z_size_w_h = local_z_size+2*halo_size,
      local_size_w_h = local_x_size_w_h * local_y_size_w_h * local_z_size_w_h,
      i_min_w_h = local_grid.I_Min_With_Halo(), 
      i_max_w_h = local_grid.I_Max_With_Halo(),
      j_min_w_h = local_grid.J_Min_With_Halo(), 
      j_max_w_h = local_grid.J_Max_With_Halo(),
      k_min_w_h = local_grid.K_Min_With_Halo(), 
      k_max_w_h = local_grid.K_Max_With_Halo(),
      global_x_size = local_x_size*mpi_driver->num_procs[0],
      global_y_size = local_y_size*mpi_driver->num_procs[1],
      global_z_size = local_z_size*mpi_driver->num_procs[2],
      ret = 1;
  int proc_coords[mpi_driver->num_dimensions];
  // 1-3 components - vector_3d, 4th - f
  T local_grid_and_F[local_x_size_w_h][local_y_size_w_h][local_z_size_w_h][4];
  MPI_Status status; 

  // create global grid and function arrays
  global_grid = new ARRAY_3D<VECTOR_3D<T> >(1, global_x_size, 1,global_y_size,
					    1, global_z_size, halo_size);
  global_F = new ARRAY_3D<T>(1, global_x_size, 1, global_y_size,
			     1, global_z_size, halo_size);
  int global_size_w_h = global_grid->Total_Size_With_Halo();
  T global_raw_array[global_size_w_h][4];

  if(mpi_driver->my_rank == 0){
    //embed local array into global for proc# = 0
    for(int i=i_min_w_h; i<=i_max_w_h; i++)
      for(int j=j_min_w_h; j<=j_max_w_h; j++)	
	for(int k=k_min_w_h; k<=k_max_w_h; k++){
	  int ig = local_x_size * mpi_driver->my_coords_in_grid[0] + i,
	      jg = local_y_size * mpi_driver->my_coords_in_grid[1] + j,
 	      kg = local_z_size * mpi_driver->my_coords_in_grid[2] + k;     
	  (*global_grid)(ig,jg,kg) = local_grid(i,j,k);
	  (*global_F)(ig,jg,kg) = local_F(i,j,k);
    }    
    //embed local array into global(located on proc#0) for proc# > 0
    for(int source = 1; source < mpi_driver->total_procs; source++){
      //receive coordinates in Cartesian grid for proc# = source
      int tag = 0;
      MPI_Recv(proc_coords, mpi_driver->num_dimensions, MPI_INT, source, tag, 
	       mpi_driver->Get_Grid_Communicator(), &status);
      //receive local array from proc# = source
      tag = 1;
      MPI_Recv(local_grid_and_F, 4*local_size_w_h, MPI_DOUBLE, source, tag, 
	       mpi_driver->Get_Grid_Communicator(), &status);
      for(int i=i_min_w_h; i<=i_max_w_h; i++)
	for(int j=j_min_w_h; j<=j_max_w_h; j++)	
	  for(int k=k_min_w_h; k<=k_max_w_h; k++){
	    int ig = local_x_size * proc_coords[0] + i,
	        jg = local_y_size * proc_coords[1] + j,
	        kg = local_z_size * proc_coords[2] + k,
	        il = i-i_min_w_h, jl = j-j_min_w_h, kl = k-k_min_w_h;
	    (*global_grid)(ig,jg,kg) = VECTOR_3D<T>(
					     local_grid_and_F[il][jl][kl][0],
					     local_grid_and_F[il][jl][kl][1],
					     local_grid_and_F[il][jl][kl][2] );
	    (*global_F)(ig,jg,kg) = local_grid_and_F[il][jl][kl][3];
	  }
    }
    //mpi_driver->Write_Local_Array_To_Disk("interpolant_grid", *global_grid);
    //mpi_driver->Write_Local_Array_To_Disk("interpolant_F", *global_F);

    // populate global flat array for broadcasting to all procs
    for(int i=0;i<global_size_w_h;i++){
      global_raw_array[i][0] = global_grid->Raw_Array(i).x;
      global_raw_array[i][1] = global_grid->Raw_Array(i).y;
      global_raw_array[i][2] = global_grid->Raw_Array(i).z;
      global_raw_array[i][3] = global_F->Raw_Array(i);
    }
  }else{//my_rank !=0
    for(int i=i_min_w_h; i<=i_max_w_h; i++)
      for(int j=j_min_w_h; j<=j_max_w_h; j++)	
	for(int k=k_min_w_h; k<=k_max_w_h; k++){ 	  
	  int il = i-i_min_w_h, jl = j-j_min_w_h, kl = k-k_min_w_h;
	  local_grid_and_F[il][jl][kl][0] = local_grid(i,j,k).x;
	  local_grid_and_F[il][jl][kl][1] = local_grid(i,j,k).y;
	  local_grid_and_F[il][jl][kl][2] = local_grid(i,j,k).z;
	  local_grid_and_F[il][jl][kl][3] = local_F(i,j,k);
	}
    int dest = 0, tag = 0;
    MPI_Send(mpi_driver->my_coords_in_grid, mpi_driver->num_dimensions, MPI_INT,
             dest, tag, mpi_driver->Get_Grid_Communicator());
    tag = 1;
    MPI_Send(local_grid_and_F, 4*local_size_w_h, MPI_DOUBLE, dest, tag, 
	     mpi_driver->Get_Grid_Communicator());
  }  

  // broadcast aggregated global grid and F arrays to all procs from 0
  MPI_Bcast(&global_raw_array, 4*global_size_w_h, MPI_DOUBLE, 0, 
	    mpi_driver->Get_Grid_Communicator());
  // copy flat array to ARRAY_3D data structures on all procs > 0
  if(mpi_driver->my_rank)
    for(int i=0;i<global_size_w_h;i++){
      global_grid->Raw_Array(i) = VECTOR_3D<T>(global_raw_array[i][0], 
                                global_raw_array[i][1], global_raw_array[i][2]);
      global_F->Raw_Array(i) = global_raw_array[i][3];
    }   
  if(mpi_driver->my_rank==1){
    mpi_driver->Write_Local_Array_To_Disk("interpolant_grid", *global_grid);
    mpi_driver->Write_Local_Array_To_Disk("interpolant_F", *global_F);
  }
  mpi_driver->Syncronize_All_Procs();
}
//*****************************************************************************
// Destructor
//*****************************************************************************
template<class T>
INTERPOLANT<T>::~INTERPOLANT() 
{  
  delete global_grid; delete global_F;
}
//*****************************************************************************
//*****************************************************************************
template<class T>
void INTERPOLANT<T>::Interpolate_On_New_Grid(
			ARRAY_3D<VECTOR_3D<T> >& new_grid, ARRAY_3D<T>& new_f)
{
  for(int i=new_f.I_Min_With_Halo(); i<=new_f.I_Max_With_Halo(); i++)
    for(int j=new_f.J_Min_With_Halo(); j<=new_f.J_Max_With_Halo(); j++)
      for(int k=new_f.K_Min_With_Halo(); k<=new_f.K_Max_With_Halo(); k++){
	VECTOR_3D<int> node1(-100,-100,-100), node2(node1), nodeNAN(node1);
	//cout<<"I="<<i<<",J="<<j<<",K="<<k<<endl;
	int kk = 1;
	Find_Confining_Cell_On_Global_Grid(kk, new_grid(i,j,k), node1, node2);
	if(node1 == nodeNAN || node2 == nodeNAN)
	  cout<<"No cell found:"<<"I="<<i<<",J="<<j<<",K="<<k<<endl;
	if((node1 != nodeNAN) && (node2 != nodeNAN))
	  if(node1 == node2) {//right no the existing node
	    //cout<<"Node interpolation, node1="<<node1<<endl;
	    new_f(i,j,k) = (*global_F)(node1.x, node1.y, node1.z);
	  }else
	    new_f(i,j,k) = Interpolate_Function_In_Quad(new_grid(i,j,k),
	                                                          node1, node2);
      }
}
//***************************************************************************** 
// Finds the cell that contains 'location' on the global_grid
//*****************************************************************************
template<class T>
void INTERPOLANT<T>::Find_Confining_Cell_On_Global_Grid(
    int k, VECTOR_3D<T>& location, VECTOR_3D<int>& node1, VECTOR_3D<int>& node2)
{
  int imin_w_h = global_grid->I_Min_With_Halo(), 
      imax_w_h = global_grid->I_Max_With_Halo(),
      jmin_w_h = global_grid->J_Min_With_Halo(), 
      jmax_w_h = global_grid->J_Max_With_Halo();
  for(int i=imin_w_h; i<=imax_w_h; i++)
    for(int j=jmin_w_h; j<=jmax_w_h; j++){      
      VECTOR_3D<T> lower_left  = (*global_grid)(i,j,k); 
      T eps = 1e-15;
      if( (fabs(lower_left.x)-eps < fabs(location.x)) && 
	  (fabs(lower_left.x)+eps > fabs(location.x)) && 
	  (fabs(lower_left.y)-eps < fabs(location.y)) && 
          (fabs(lower_left.y)+eps > fabs(location.y)) )
	{cout<<"Node"<<endl; node1 = node2 = VECTOR_3D<int>(i,j,k); return;}
      
      // for vertical moving grid only
      if(j<jmax_w_h){
	T upper_right_y = (*global_grid)(i,j+1,k).y;
        if(fmin(lower_left.y, upper_right_y) <= location.y && 
	   fmax(lower_left.y, upper_right_y) >= location.y) {
	  //cout<<"Cell"<<endl;
	  node1 = VECTOR_3D<int>(i,j,k); node2 = VECTOR_3D<int>(i,j+1,k);
	  return;
	}
      /*      
	if(i<imax_w_h && j<jmax_w_h){
	VECTOR_3D<T> upper_right = (*global_grid)(i+1,j+1,k);
        if((fmin(lower_left.x, upper_right.x) <= location.x && 
	    fmax(lower_left.x, upper_right.x) >= location.x   ) &&
	   (fmin(lower_left.y, upper_right.y) <= location.y && 
	    fmax(lower_left.y, upper_right.y) >= location.y)  ) {
	  //cout<<"Cell"<<endl;
	  node1 = VECTOR_3D<int>(i,j,k); node2 = VECTOR_3D<int>(i+1,j+1,k);
	  return;
	}
      */
      }
    }
}
//*****************************************************************************
// Finds the value of function in location P (belonging to the quad:node1,node2)
//*****************************************************************************
template<class T>
T INTERPOLANT<T>::Interpolate_Function_In_Quad(VECTOR_3D<T>& P,
                                   VECTOR_3D<int>& node1, VECTOR_3D<int>& node2)
{
  // for vertical moving (in y) grid only
  assert(node1.x == node2.x && node1.z == node2.z);
  T Ay = (*global_grid)(node1.x, node1.y, node1.z).y, 
    By = (*global_grid)(node2.x, node2.y, node2.z).y,
    FA = (*global_F)(node1.x, node1.y, node1.z), 
    FB = (*global_F)(node2.x, node2.y, node2.z),
    AP = fabs(Ay-P.y), AB = fabs(Ay-By);
    assert(AB);
  T wA = (AB-AP)/AB;
  assert(wA>=0. && wA<=1.);
  return FA*wA + FB*(1.-wA);
  /*
  // Interpolating in one of two Triangles, forming the Quad
  VECTOR_3D<T> A = (*global_grid)(node1.x, node1.y, node1.z), 
               B = (*global_grid)(node2.x, node1.y, node1.z),
               C = (*global_grid)(node2.x, node2.y, node1.z),
               D = (*global_grid)(node1.x, node2.y, node1.z);

  T FA = (*global_F)(node1.x, node1.y, node1.z), 
    FB = (*global_F)(node2.x, node1.y, node1.z),
    FC = (*global_F)(node2.x, node2.y, node1.z),
    FD = (*global_F)(node1.x, node2.y, node1.z);  

  T AP = (A-P).Magnitude(), BP = (B-P).Magnitude(), 
    CP = (C-P).Magnitude(), DP = (D-P).Magnitude();
  if(AP>BP && AP>CP && AP>DP) 
    return Interpolate_Function_In_Triangle(P, B,C,D, FB,FC,FD);
  else if(BP>AP && BP>CP && BP>DP) 
    return Interpolate_Function_In_Triangle(P, A,C,D, FA,FC,FD);
  else if(CP>AP && CP>BP && CP>DP) 
    return Interpolate_Function_In_Triangle(P, A,B,D, FA,FB,FD);
  else if(DP>AP && DP>BP && DP>CP) 
    return Interpolate_Function_In_Triangle(P, A,B,C, FA,FB,FC);
  else
    assert(false);
  */
  /*
  // Interpolating in Quad
  T aq = VECTOR_3D<T>::Cross_Product(A-D, C-B).Magnitude(),
    bq = VECTOR_3D<T>::Cross_Product(P-A, A-D+C-B).Magnitude() + 
         VECTOR_3D<T>::Cross_Product(A-D, B-A).Magnitude(),
    cq = VECTOR_3D<T>::Cross_Product(P-A, B-A).Magnitude(),
    t1 = (-bq + sqrt(bq*bq - 4.*aq*cq))/2.*aq,
    t2 = (-bq - sqrt(bq*bq - 4.*aq*cq))/2.*aq,
    t = ((t1>=0 && t1<=1) ? t1 : t2),
    s = (P - A + t * (A - D)).x / (B - A + t * (A - D + C - B)).x;
  //if(!(t1>=0 && t1<=1) && !(t2>=0 && t2<=1))
  //  cout <<"T1="<<t1<<", T2="<<t2<<endl;
  assert( (t1>=0 && t1<=1) || (t2>=0 && t2<=1) );

  T F_A = (*global_F)(node1.x, node1.y, node1.z), 
    F_B = (*global_F)(node2.x, node1.y, node1.z),
    F_D = (*global_F)(node1.x, node2.y, node1.z),
    F_C = (*global_F)(node2.x, node2.y, node1.z),
    F_P = F_A + (F_B-F_A)*s,
    F_Q = F_D + (F_C-F_D)*s,
    F = F_P + (F_Q-F_P)*t;
  //cout<<"F = "<<F<<",t="<<t<<",s="<<s<<endl;
  return F;*/
}
//*****************************************************************************
template<class T>
T INTERPOLANT<T>::Interpolate_Function_In_Triangle(VECTOR_3D<T>& p,
         VECTOR_3D<T>& v1, VECTOR_3D<T>& v2, VECTOR_3D<T>& v3, T f1, T f2, T f3)
{
  VECTOR_3D<T> u=v2-v1, v=v3-v1, w=p-v1;
    T u_u=VECTOR_3D<T>::Dot_Product(u,u), v_v=VECTOR_3D<T>::Dot_Product(v,v),
      u_v=VECTOR_3D<T>::Dot_Product(u,v), u_w=VECTOR_3D<T>::Dot_Product(u,w),
      v_w=VECTOR_3D<T>::Dot_Product(v,w);
    T denom = u_u*v_v - u_v*u_v, one_over_denom;
    if(fabs(denom)>1e-16) one_over_denom=1./denom;else one_over_denom=(T)1e16;
    T a=(v_v*u_w-u_v*v_w)*one_over_denom, b=(u_u*v_w-u_v*u_w)*one_over_denom;
    return f1*(1-a-b) + f2*a + f3*b;
}
//*****************************************************************************
#endif
