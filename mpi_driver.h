//*****************************************************************************
// Note: Class for MPI framework management
//*****************************************************************************
#ifndef __MPI_DRIVER__
#define __MPI_DRIVER__

#include <mpi.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "array_3d.h"
#include "vector_3d.h"
#include "parameters.h"
#include <sys/stat.h>

using namespace std;

template<class T=double>
class MPI_DRIVER
{
  public:
    MPI_DRIVER(PARAMETERS<T>& p) 
      : output_dir(p.output_dir)
    {Initialize(p.argc,p.argv, p.num_cpu_x, p.num_cpu_y, p.num_cpu_z, 
        p.periodic_in_x, p.periodic_in_y, p.periodic_in_z, 
        p.num_local_nodes_x, p.num_local_nodes_y, p.num_local_nodes_z);}

    ~MPI_DRIVER()
    {MPI_Finalize();
      delete [] periodic; delete [] num_procs; delete [] my_coords_in_grid;
      delete [] local_grid_lower_bound; delete [] local_grid_upper_bound;}

    void Initialize(int argc, char* argv[], 
        int num_procs_in_x, int num_procs_in_y, int num_procs_in_z, 
        bool x_periodic, bool y_periodic, bool z_periodic, 
        int local_size_x, int local_size_y, int local_size_z);
    void Initialize_MPI(int argc, char* argv[], int num_x, int num_y, int num_z);
    void Exchange_Ghost_Values_For_Scalar_Field(ARRAY_3D<T>& scalar, 
        int Halo_size=-1);
    void Exchange_Ghost_Values_For_Vector_Field(ARRAY_3D<VECTOR_3D<T> >& vector);
    void Replace_With_Sum_On_All_Procs(T& value);
    void Replace_With_Max_Value_Among_All_Procs(T& value);
    void Replace_With_Max_Value_Among_All_Procs(int& value);

    int Write_Global_Array_To_Disk(string a_name, ARRAY_3D<T>& a, 
        int timestep=0, bool save_halo = false);
    //int Write_Global_Array_To_Disk(string a_name, ARRAY_2D<T>& a);
    int Write_Global_Array_To_Disk(string va_name, ARRAY_3D<VECTOR_3D<T> >& va, 
        int timestep = 0, bool save_halo = false);

    int Read_Global_Array_From_Disk(string filename, ARRAY_3D<VECTOR_3D<T> >& a);

    //int Write_Global_Array_To_Disk(string va_name, ARRAY_2D<VECTOR_3D<T> >& va);
    int Write_Local_Array_To_Disk(string a_name, ARRAY_3D<T>& a, 
        bool add_proc_number = true, int timestep = 0);
    int Write_Local_Array_To_Disk(string a_name, ARRAY_2D<T>& a, int timestep = 0)
    {ARRAY_3D<T> a3d(a); 
      return Write_Local_Array_To_Disk(a_name, a3d,true,timestep);}
    int Write_Local_Array_To_Disk(string va_name, ARRAY_3D<VECTOR_3D<T> >& va, 
        int timestep = 0);
    int Write_Local_Array_To_Disk(string va_name, ARRAY_2D<VECTOR_3D<T> >& va, 
        int timestep = 0)
    {ARRAY_3D<VECTOR_3D<T> > va3d(va); 
      return Write_Local_Array_To_Disk(va_name, va3d, timestep);}
    int Write_Local_Array_To_Disk(string a_name, T* a, int size, int timestep=0);
    // binary write funcs
    int Write_Binary_Local_Array(string a_name, ARRAY_3D<T>& a);
    void Write_Binary_Local_Array(ofstream& output,ARRAY_3D<T>& a);  
    void Write_Binary_Local_Array(ofstream& output,ARRAY_3D<VECTOR_3D<T> >& va);  
    int Read_Binary_Local_Array(string a_name, ARRAY_3D<T>& a);
    void Read_Binary_Local_Array(ifstream& input, ARRAY_3D<T>& a); 
    void Read_Binary_Local_Array(ifstream& input, ARRAY_3D<VECTOR_3D<T> >& va);
    void Syncronize_All_Procs() {MPI_Barrier(grid_comm);}
    MPI_Comm& Get_Grid_Communicator() {return grid_comm;}


    int my_rank, total_procs, num_dimensions;
    int *num_procs, *my_coords_in_grid;
    int west_proc, east_proc, suth_proc, nrth_proc, back_proc, frnt_proc;
    int *local_grid_lower_bound, *local_grid_upper_bound;

  private:
    MPI_Comm grid_comm;
    MPI_Status status; 
    int *periodic;
    string output_dir, proc_coords_filename;

    int Save_Procs_Coordinates_In_File(string& output_filename);
    void Create_Output_Directory_If_Missing();
};

//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Initialize(int argc, char* argv[], int num_procs_in_x, int num_procs_in_y, int num_procs_in_z, bool x_periodic, bool y_periodic, bool z_periodic, int local_size_x, int local_size_y, int local_size_z)
{
  num_dimensions = 3;
  periodic    = new int[num_dimensions];
  periodic[0] = (x_periodic ? 1 : 0); 
  periodic[1] = (y_periodic ? 1 : 0); 
  periodic[2] = (z_periodic ? 1 : 0);

  my_coords_in_grid = new int[num_dimensions];
  for(int i=0;i<num_dimensions;i++) my_coords_in_grid[i] = 0;

  num_procs  = new int[num_dimensions];
  num_procs[0] = num_procs_in_x; 
  num_procs[1] = num_procs_in_y; 
  num_procs[2] = num_procs_in_z;
  total_procs = num_procs_in_x * num_procs_in_y * num_procs_in_z;
  // local grid index bounds
  local_grid_lower_bound = new int[num_dimensions]; 
  local_grid_upper_bound = new int[num_dimensions];

  proc_coords_filename = output_dir + "cpu_coords.txt";
  Initialize_MPI(argc,argv, local_size_x, local_size_y, local_size_z);
}
//*****************************************************************************
// Sets up a grid communicator and saves coords for each proc to a file
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Initialize_MPI(int argc, char* argv[], 
    int local_size_x, int local_size_y, int local_size_z)
{
  int reorder = 1, total_procs_mpi;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total_procs_mpi);
  if(total_procs != total_procs_mpi) 
    cout<<"WARNING: number of procs does not match (required: "<<total_procs
      <<", available:"<<total_procs_mpi<<")"<<endl;

  MPI_Cart_create(MPI_COMM_WORLD, num_dimensions, num_procs, periodic, 
      reorder,&grid_comm);
  MPI_Cart_coords(grid_comm, my_rank, num_dimensions, my_coords_in_grid);

  // find ranks of neighbors
  MPI_Cart_shift(grid_comm, 0, 1, &west_proc, &east_proc);
  MPI_Cart_shift(grid_comm, 1, 1, &suth_proc, &nrth_proc);  
  MPI_Cart_shift(grid_comm, 2, 1, &back_proc, &frnt_proc);
  // save lower and upper index bounds for local domains
  if(west_proc!=MPI_PROC_NULL) local_grid_lower_bound[0]=0; 
  else local_grid_lower_bound[0]=1;
  if(east_proc!=MPI_PROC_NULL) local_grid_upper_bound[0]=local_size_x; 
  else local_grid_upper_bound[0]=local_size_x-1;
  if(suth_proc!=MPI_PROC_NULL) local_grid_lower_bound[1]=0; 
  else local_grid_lower_bound[1]=1;
  if(nrth_proc!=MPI_PROC_NULL) local_grid_upper_bound[1]=local_size_y; 
  else local_grid_upper_bound[1]=local_size_y-1;
  if(back_proc!=MPI_PROC_NULL) local_grid_lower_bound[2]=0; 
  else local_grid_lower_bound[2]=1;
  if(frnt_proc!=MPI_PROC_NULL) local_grid_upper_bound[2]=local_size_z; 
  else local_grid_upper_bound[2]=local_size_z-1;
  // save processor Cartesian grid cooordinates in a file
  Save_Procs_Coordinates_In_File(proc_coords_filename);
}
//*****************************************************************************
// Saving grid communicator coordinates in a text file 'output_filename'
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Save_Procs_Coordinates_In_File(string& output_filename)
{
  int proc_value, proc_coords[num_dimensions];

  if(my_rank==0){
    Create_Output_Directory_If_Missing(); 
    ofstream output_file(output_filename.c_str(), ios::out);
    if(!output_file) {
      cout<<"ERROR: could not open file for writing"<<endl;
      return 0;
    }
    // table header
    output_file.width(7);
    output_file <<right<<"CPU_NUM";
    output_file.width(14);
    output_file <<right<<"CPU_X_COORD";
    output_file.width(14);
    output_file <<right<<"CPU_Y_COORD";
    output_file.width(14);
    output_file <<right<<"CPU_Z_COORD" << endl;
    // output coords for Proc# 0
    output_file.width(7);
    output_file <<right<< my_rank;
    for(int i=0; i<num_dimensions; i++){
      output_file.width(14); 
      output_file <<right<< my_coords_in_grid[i];
    }
    output_file << endl;  
    //cout<<"MPI_PROC_NULL="<<MPI_PROC_NULL<<endl;      
    //cout<<"Proc#"<<my_rank<<",E.P="<<east_proc<<",W.P="<<west_proc
    //	  <<",S.P="<<suth_proc<<",N.P="<<nrth_proc 
    //    <<",B.P="<<back_proc<<",F.P="<<frnt_proc <<endl;
    // output coords for Proc# 1..total_procs
    for(int source = 1; source < total_procs; source++){
      int tag = 0;
      MPI_Recv(&proc_value, 1, MPI_INT, source, tag, grid_comm, &status);
      tag = 1;
      MPI_Recv(proc_coords, num_dimensions, MPI_INT, source, tag, 
          grid_comm, &status);
      output_file.width(7);     
      output_file <<right<< proc_value;
      for(int i=0; i<num_dimensions; i++){
        output_file.width(14);
        output_file <<right<< proc_coords[i];
      }
      output_file << endl;

      int ep,wp,sp,np,bp,fp;
      MPI_Recv(&ep,1,MPI_INT,source,2,grid_comm, &status);
      MPI_Recv(&wp,1,MPI_INT,source,3,grid_comm, &status);
      MPI_Recv(&sp,1,MPI_INT,source,4,grid_comm, &status);
      MPI_Recv(&np,1,MPI_INT,source,5,grid_comm, &status);
      MPI_Recv(&bp,1,MPI_INT,source,6,grid_comm, &status);
      MPI_Recv(&fp,1,MPI_INT,source,7,grid_comm, &status);
      //cout<<"Proc#"<<proc_value <<",E.P="<<ep<<",W.P="<<wp
      //  <<",S.P="<<sp<<",N.P="<<np <<",B.P="<<bp<<",F.P="<<fp <<endl;
    }//for: source 
    output_file.close();
  }//if: my_rank==0
  else{// send my_rank and my_coords_in_grid[] to Proc# 0
    int dest = 0;
    int tag = 0;
    MPI_Send(&my_rank, 1, MPI_INT, dest, tag, grid_comm);
    tag = 1;
    MPI_Send(my_coords_in_grid,num_dimensions,MPI_INT,dest,tag,grid_comm);
    MPI_Send(&east_proc,1,MPI_INT,dest,2,grid_comm);
    MPI_Send(&west_proc,1,MPI_INT,dest,3,grid_comm);
    MPI_Send(&suth_proc,1,MPI_INT,dest,4,grid_comm);
    MPI_Send(&nrth_proc,1,MPI_INT,dest,5,grid_comm);
    MPI_Send(&back_proc,1,MPI_INT,dest,6,grid_comm);
    MPI_Send(&frnt_proc,1,MPI_INT,dest,7,grid_comm);
  }
  return 1;
}
//*****************************************************************************
// Outputs to disk global ARRAY_3D<VECTOR_3D<T> > (w/o halo) on proc#0
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Write_Global_Array_To_Disk(string va_name, 
    ARRAY_3D<VECTOR_3D<T> >& va, int timestep, bool save_halo)
{
  int ret = 1;
  ARRAY_3D<T> vx(va, 1); ARRAY_3D<T> vy(va, 2); ARRAY_3D<T> vz(va, 3);
  ret = min(ret,Write_Global_Array_To_Disk(va_name+"_x",vx,timestep,save_halo));
  ret = min(ret,Write_Global_Array_To_Disk(va_name+"_y",vy,timestep,save_halo));
  ret = min(ret,Write_Global_Array_To_Disk(va_name+"_z",vz,timestep,save_halo));
  return ret; // ret=0 in case either one of three writes fails
} 
//*****************************************************************************
// Outputs to disk global ARRAY_3D<T> on proc#0
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Write_Global_Array_To_Disk(string a_name, ARRAY_3D<T>& a, 
    int timestep, bool save_halo)
{
  int local_x_size_nh = a.I_Size(), local_y_size_nh = a.J_Size(),
      local_z_size_nh = a.K_Size(),
      h = (save_halo ? a.Halo_Size() : 0), hh = 2*h,
      local_x_size = local_x_size_nh + hh, local_y_size = local_y_size_nh + hh, 
      local_z_size = local_z_size_nh + hh,
      local_size = local_x_size * local_y_size * local_z_size,
      glob_x_size  = local_x_size_nh*num_procs[0] + hh,
      glob_y_size  = local_y_size_nh*num_procs[1] + hh,
      glob_z_size  = local_z_size_nh*num_procs[2] + hh,
      a_i_min = (save_halo ? a.I_Min_With_Halo() : a.I_Min()),
      a_j_min = (save_halo ? a.J_Min_With_Halo() : a.J_Min()),
      a_k_min = (save_halo ? a.K_Min_With_Halo() : a.K_Min()),
      a_i_max = (save_halo ? a.I_Max_With_Halo() : a.I_Max()),
      a_j_max = (save_halo ? a.J_Max_With_Halo() : a.J_Max()),
      a_k_max = (save_halo ? a.K_Max_With_Halo() : a.K_Max()),
      ret = 1,
      proc_coords[num_dimensions];
  T local_array[local_x_size][local_y_size][local_z_size];

  if(my_rank==0){
    ARRAY_3D<T> global_a(1,glob_x_size, 1,glob_y_size, 1,glob_z_size, 0);//no h
    //embed local array into global for proc# = 0
    for(int i=a_i_min; i<=a_i_max; i++)
      for(int j=a_j_min; j<=a_j_max; j++)	
        for(int k=a_k_min; k<=a_k_max; k++){
          int i_g = h + local_x_size_nh*my_coords_in_grid[0] + i,
              j_g = h + local_y_size_nh*my_coords_in_grid[1] + j,    
              k_g = h + local_z_size_nh*my_coords_in_grid[2] + k;
          global_a(i_g,j_g,k_g) = a(i,j,k);
        }
    //embed local array into global(located on proc#0) for proc# > 0
    for(int source = 1; source < total_procs; source++){
      //receive coordinates in Cartesian grid for proc# = source
      int tag = 0;
      MPI_Recv(proc_coords, num_dimensions, MPI_INT, source, tag, 
          grid_comm, &status);
      //receive local array from proc# = source
      tag = 1;
      MPI_Recv(local_array, local_size, MPI_DOUBLE, source, tag, 
          grid_comm, &status);
      for(int i=a_i_min; i<=a_i_max; i++)
        for(int j=a_j_min; j<=a_j_max; j++)	
          for(int k=a_k_min; k<=a_k_max; k++){
            int i_g = h + local_x_size_nh*proc_coords[0] + i,
                j_g = h + local_y_size_nh*proc_coords[1] + j,
                k_g = h + local_z_size_nh*proc_coords[2] + k;
            global_a(i_g,j_g,k_g)= local_array[i-a_i_min][j-a_j_min][k-a_k_min];
          }
    }
    //write global array on proc#0
    ret = Write_Local_Array_To_Disk(a_name,global_a,false,timestep);
  }else{//my_rank !=0
    for(int i=a_i_min; i<=a_i_max; i++)
      for(int j=a_j_min; j<=a_j_max; j++)	
        for(int k=a_k_min; k<=a_k_max; k++)
          local_array[i-a_i_min][j-a_j_min][k-a_k_min] = a(i,j,k);   
    int dest = 0, tag = 0;
    MPI_Send(my_coords_in_grid,num_dimensions,MPI_INT,dest,tag,grid_comm);
    tag = 1;
    MPI_Send(local_array,local_size,MPI_DOUBLE,dest,tag,grid_comm);
  }
  return ret;
}
//*****************************************************************************
// Reads global ARRAY_3D<VECTOR_3D<T> > (w/ halo) on proc#0 and distributes it
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Read_Global_Array_From_Disk(string filename, 
    ARRAY_3D<VECTOR_3D<T> >& a)
{
  int loc_i_size = a.I_Max() - a.I_Min() + 1, 
      loc_j_size = a.J_Max() - a.J_Min() + 1,
      loc_k_size = a.K_Max() - a.K_Min() + 1, halo = a.Halo_Size(),
      proc_coords[num_dimensions],
      loc_array_size = 
        3*(loc_i_size+2*halo) * (loc_j_size+2*halo) * (loc_k_size+2*halo);
  T loc_array[loc_i_size+2*halo][loc_j_size+2*halo][loc_k_size+2*halo][3];

  if(my_rank==0){
    int glob_i_size = loc_i_size*num_procs[0],
        glob_j_size = loc_j_size*num_procs[1],
        glob_k_size = loc_k_size*num_procs[2];
    ARRAY_3D<VECTOR_3D<T> > glob_a(1, glob_i_size, 1, glob_j_size, 
        1, glob_k_size, halo);

    // Read global array from file
    ifstream infile(filename.c_str(), ios::in|ios::binary);
    if(infile.is_open()){
      int size_i, size_j, size_k;
      infile.read((char*)(&size_i), sizeof(int));
      assert(size_i == glob_a.I_Size_With_Halo());
      infile.read((char*)(&size_j), sizeof(int)); 
      assert(size_j == glob_a.J_Size_With_Halo());
      infile.read((char*)(&size_k), sizeof(int)); 
      assert(size_k == glob_a.K_Size_With_Halo());

      for(int i=glob_a.I_Min_With_Halo(); i<=glob_a.I_Max_With_Halo(); i++)
        for(int j=glob_a.J_Min_With_Halo(); j<=glob_a.J_Max_With_Halo(); j++)
          for(int k=glob_a.K_Min_With_Halo(); k<=glob_a.K_Max_With_Halo(); k++){
            infile.read((char*)(&glob_a(i,j,k).x),sizeof(double));
            infile.read((char*)(&glob_a(i,j,k).y),sizeof(double));
            infile.read((char*)(&glob_a(i,j,k).z),sizeof(double));
          }
      infile.close();
    }else {cout<<"ERROR: Cannot open the global array file"<<endl; return 0;}

    //extract local array from global for proc# = 0
    for(int i=a.I_Min_With_Halo(); i<=a.I_Max_With_Halo(); i++)
      for(int j=a.J_Min_With_Halo(); j<=a.J_Max_With_Halo(); j++)	
        for(int k=a.K_Min_With_Halo(); k<=a.K_Max_With_Halo(); k++){
          int i_g = loc_i_size * my_coords_in_grid[0]+i,//index into glob arr
              j_g = loc_j_size * my_coords_in_grid[1]+j,
              k_g = loc_k_size * my_coords_in_grid[2]+k;     
          a(i,j,k) = glob_a(i_g,j_g,k_g);
        }    
    //extract&send local array from global(located on proc#0) for proc# > 0
    for(int source = 1; source < total_procs; source++){
      //receive coordinates in Cartesian grid for proc# = source
      int tag = 0;
      MPI_Recv(proc_coords, num_dimensions, MPI_INT, source, tag, 
          grid_comm, &status);
      //send local_array to proc# = source
      for(int i=a.I_Min_With_Halo(); i<=a.I_Max_With_Halo(); i++)
        for(int j=a.J_Min_With_Halo(); j<=a.J_Max_With_Halo(); j++)	
          for(int k=a.K_Min_With_Halo(); k<=a.K_Max_With_Halo(); k++){
            int i_g = loc_i_size * proc_coords[0]+i, //index into global array
                j_g = loc_j_size * proc_coords[1]+j,
                k_g = loc_k_size * proc_coords[2]+k,     
                i_l = i-a.I_Min_With_Halo(), 
                j_l = j-a.J_Min_With_Halo(), 
                k_l = k-a.K_Min_With_Halo();
            loc_array[i_l][j_l][k_l][0] = glob_a(i_g,j_g,k_g).x;
            loc_array[i_l][j_l][k_l][1] = glob_a(i_g,j_g,k_g).y;
            loc_array[i_l][j_l][k_l][2] = glob_a(i_g,j_g,k_g).z;
          }   
      tag = 1;    
      MPI_Send(loc_array, loc_array_size, MPI_DOUBLE, source, tag, grid_comm);
    }
  }else{//my_rank !=0
    //send coordinates in the global array
    int dest = 0, tag = 0;
    MPI_Send(my_coords_in_grid, num_dimensions, MPI_INT, dest, tag, grid_comm);
    //receive local_array from proc#0
    tag = 1;
    MPI_Recv(loc_array, loc_array_size, MPI_DOUBLE,dest,tag,grid_comm,&status);
    //populate the local array
    for(int i=a.I_Min_With_Halo(); i<=a.I_Max_With_Halo(); i++)
      for(int j=a.J_Min_With_Halo(); j<=a.J_Max_With_Halo(); j++)	
        for(int k=a.K_Min_With_Halo(); k<=a.K_Max_With_Halo(); k++){
          int i_l = i-a.I_Min_With_Halo(), 
              j_l = j-a.J_Min_With_Halo(), 
              k_l = k-a.K_Min_With_Halo();
          a(i,j,k).x = loc_array[i_l][j_l][k_l][0];
          a(i,j,k).y = loc_array[i_l][j_l][k_l][1];
          a(i,j,k).z = loc_array[i_l][j_l][k_l][2];
        }
  }
  return 1;
}
/*
//*****************************************************************************
// Outputs to disk global ARRAY_2D<VECTOR_3D<T> > (w/o halo) on proc#0
//*****************************************************************************
template <class T>
int MPI_DRIVER<T>::Write_Global_Array_To_Disk(string va_name, ARRAY_2D<VECTOR_3D<T> >& va)
{
int ret = 1;
ARRAY_2D<T> vx(va, 1); ARRAY_2D<T> vy(va, 2); ARRAY_2D<T> vz(va, 3);
ret = min(ret, Write_Global_Array_To_Disk(va_name+"_x", vx));
ret = min(ret, Write_Global_Array_To_Disk(va_name+"_y", vy));
ret = min(ret, Write_Global_Array_To_Disk(va_name+"_z", vz));
return ret; // ret=0 in case either one of three writes fails
}
//*****************************************************************************
// Outputs to disk global ARRAY_2D<T> (w/o halo) on proc#0
//*****************************************************************************
template <class T>
int MPI_DRIVER<T>::Write_Global_IJ_Array_To_Disk(string a_name, ARRAY_2D<T>& a)
{
if(proc_coords[2] == 0){ // z coord of processor in Cartesian topology = 0
int local_x_size = a.I_Max() - a.I_Min() + 1, 
local_y_size = a.J_Max() - a.J_Min() + 1,
local_size = local_x_size * local_y_size;
int ret=1;
int proc_coords[num_dimensions];
T local_array[local_x_size][local_y_size];

if(my_rank==0){
ARRAY_2D<T> global_a(1, local_x_size*num_procs[0], 
1, local_y_size*num_procs[1], 0); //no halo
//embed local array into global for proc# = 0
for(int i=a.I_Min();i<=a.I_Max();i++)
for(int j=a.J_Min();j<=a.J_Max();j++){
int ii = local_x_size * my_coords_in_grid[0]+i,//ind into global array
jj = local_y_size * my_coords_in_grid[1]+j,    
global_a(ii,jj) = a(i,j);
}    
//embed local array into global(located on proc#0) for proc# > 0
for(int source = 1; source < total_procs; source++){
//receive coordinates in Cartesian grid for proc# = source
int tag = 0;
MPI_Recv(proc_coords, num_dimensions, MPI_INT, source, tag, 
grid_comm, &status);
//receive local array from proc# = source
tag = 1;
MPI_Recv(local_array, local_size, MPI_DOUBLE, source, tag, 
grid_comm, &status);
for(int i=a.I_Min();i<=a.I_Max();i++)
for(int j=a.J_Min();j<=a.J_Max();j++)	
for(int k=a.K_Min();k<=a.K_Max();k++){
int ii = local_x_size * proc_coords[0]+i, //index into global array
jj = local_y_size * proc_coords[1]+j,
kk = local_z_size * proc_coords[2]+k;     
global_a(ii,jj,kk) = 
local_array[i-a.I_Min()][j-a.J_Min()][k-a.K_Min()];
}
}
//write global array on proc#0
ret = Write_Local_Array_To_Disk(a_name,global_a,false);
}else{//my_rank !=0
for(int i=a.I_Min();i<=a.I_Max();i++)
for(int j=a.J_Min();j<=a.J_Max();j++)	
for(int k=a.K_Min();k<=a.K_Max();k++)	  
local_array[i-a.I_Min()][j-a.J_Min()][k-a.K_Min()] = a(i,j,k);
int dest = 0, tag = 0;
MPI_Send(my_coords_in_grid,num_dimensions,MPI_INT,dest,tag,grid_comm);
tag = 1;
MPI_Send(local_array,local_size,MPI_DOUBLE,dest,tag,grid_comm);
}
return ret;
}
}
*/
//*****************************************************************************
// Outputs local portion of global array to a file '$a_name.#proc' on disk
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Write_Local_Array_To_Disk(string va_name, ARRAY_3D<VECTOR_3D<T> >& va, int timestep)
{
  int ret = 1;
  ARRAY_3D<T> vx(va, 1); ARRAY_3D<T> vy(va, 2); ARRAY_3D<T> vz(va, 3);
  ret = min(ret, Write_Local_Array_To_Disk(va_name+"_x", vx, true, timestep));
  ret = min(ret, Write_Local_Array_To_Disk(va_name+"_y", vy, true, timestep));
  ret = min(ret, Write_Local_Array_To_Disk(va_name+"_z", vz, true, timestep));
  return ret; // ret=0 in case either one of three writes fails
}
//*****************************************************************************
// Outputs local portion of global array to a file '$a_name.#proc' on disk
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Write_Local_Array_To_Disk(string a_name, ARRAY_3D<T>& a, bool add_proc_number, int timestep)
{
  std::stringstream filename;
  filename << output_dir;
  filename << a_name;
  if(timestep)
    filename << "_t" << timestep;
  if(add_proc_number) 
    filename << "_p" << my_rank; else filename << "_global";
  filename << ".m";

  ofstream output_file(filename.str().c_str(), ios::out);
  if(!output_file){
    cout<<"ERROR: could not open file for writing"<<endl;
    return 0;
  }
  output_file << a_name;
  if(timestep) output_file <<"_"<<timestep;
  output_file << "(:,:,:";
  if(add_proc_number) output_file <<","<<my_rank+1;
  output_file << ") = cat(3, ";

  for(int k = a.K_Min_With_Halo(); k <= a.K_Max_With_Halo(); k++){
    output_file<<"[";
    for(int j = a.J_Min_With_Halo(); j <= a.J_Max_With_Halo(); j++){
      for(int i = a.I_Min_With_Halo(); i <= a.I_Max_With_Halo(); i++){
        output_file.precision(15);
        output_file  << a(i,j,k); 
        if(i != a.I_Max_With_Halo()) output_file<<" ";
      }
      if(j != a.J_Max_With_Halo()) output_file<<"; "; else output_file<<"]"; 
      //output_file<<endl;
    }
    if(k != a.K_Max_With_Halo()) output_file<<","; else output_file<<");";
    //output_file<<endl;
  }   
  output_file.close();
  return 1;
}  
//*****************************************************************************
// Outputs local portion of global SCALAR array to a file '$a_name.#proc'
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Write_Binary_Local_Array(string a_name, ARRAY_3D<T>& a)
{
  std::stringstream filename;
  filename << output_dir << a_name << "." << my_rank; 

  ofstream output(filename.str().c_str(), ios::out | ios::binary);
  if(!output){
    cout<<"ERROR: could not open file for writing"<<endl;
    return 0;
  }
  Write_Binary_Local_Array(output, a);
  output.close();
  return 1;
}
//*****************************************************************************
// Helper: Outputs local portion of global SCALAR array to the ofstream
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Write_Binary_Local_Array(ofstream& output, ARRAY_3D<T>& a)
{
  int size=a.Total_Size_With_Halo(), imin=a.I_Min(), imax=a.I_Max(),
      jmin=a.J_Min(), jmax=a.J_Max(),kmin=a.K_Min(), kmax=a.K_Max(), 
      halo_size=a.Halo_Size();
  output.write(reinterpret_cast<char *>(&size),sizeof(int));
  output.write(reinterpret_cast<char *>(&imin),sizeof(int));
  output.write(reinterpret_cast<char *>(&imax),sizeof(int));
  output.write(reinterpret_cast<char *>(&jmin),sizeof(int));
  output.write(reinterpret_cast<char *>(&jmax),sizeof(int));
  output.write(reinterpret_cast<char *>(&kmin),sizeof(int));
  output.write(reinterpret_cast<char *>(&kmax),sizeof(int));
  output.write(reinterpret_cast<char *>(&halo_size),sizeof(int));
  output.write(reinterpret_cast<char *>(a.Raw_Array_Pointer()),sizeof(T)*size); 
} 
//*****************************************************************************
// Helper: Outputs local portion of global VECTOR array to the ofstream
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Write_Binary_Local_Array(ofstream& output, 
    ARRAY_3D<VECTOR_3D<T> >& va)
{
  ARRAY_3D<T> vx(va, 1), vy(va, 2), vz(va, 3);
  int size=va.Total_Size_With_Halo(), imin=va.I_Min(), imax=va.I_Max(),
      jmin=va.J_Min(), jmax=va.J_Max(),kmin=va.K_Min(), kmax=va.K_Max(), 
      halo_size=va.Halo_Size();
  output.write(reinterpret_cast<char *>(&size),sizeof(int));
  output.write(reinterpret_cast<char *>(&imin),sizeof(int));
  output.write(reinterpret_cast<char *>(&imax),sizeof(int));
  output.write(reinterpret_cast<char *>(&jmin),sizeof(int));
  output.write(reinterpret_cast<char *>(&jmax),sizeof(int));
  output.write(reinterpret_cast<char *>(&kmin),sizeof(int));
  output.write(reinterpret_cast<char *>(&kmax),sizeof(int));
  output.write(reinterpret_cast<char *>(&halo_size),sizeof(int));
  output.write(reinterpret_cast<char *>(vx.Raw_Array_Pointer()),sizeof(T)*size);
  output.write(reinterpret_cast<char *>(vy.Raw_Array_Pointer()),sizeof(T)*size);
  output.write(reinterpret_cast<char *>(vz.Raw_Array_Pointer()),sizeof(T)*size);
} 
//*****************************************************************************
// Loads local portion of global array from a file '$a_name.#proc' on disk
//*****************************************************************************
  template <class T>
int MPI_DRIVER<T>::Read_Binary_Local_Array(string a_name, ARRAY_3D<T>& a)
{
  std::stringstream filename;
  filename << output_dir << a_name << "." << my_rank; 

  ifstream input(filename.str().c_str(), ios::in | ios::binary);
  if(!input){
    cout<<"ERROR: could not open file for reading"<<endl;
    return 0;
  }
  Read_Binary_Local_Array(input, a);
  input.close();
  return 1;
}
//*****************************************************************************
// Helper: Outputs local portion of global VECTOR array to the ofstream
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Read_Binary_Local_Array(ifstream& input, ARRAY_3D<T>& a)
{
  int size, imin,imax,jmin,jmax,kmin,kmax,halo_size;
  input.read(reinterpret_cast<char *>(&size),sizeof(int));
  input.read(reinterpret_cast<char *>(&imin),sizeof(int));
  input.read(reinterpret_cast<char *>(&imax),sizeof(int));
  input.read(reinterpret_cast<char *>(&jmin),sizeof(int));
  input.read(reinterpret_cast<char *>(&jmax),sizeof(int));
  input.read(reinterpret_cast<char *>(&kmin),sizeof(int));
  input.read(reinterpret_cast<char *>(&kmax),sizeof(int));
  input.read(reinterpret_cast<char *>(&halo_size),sizeof(int));
  a.Delete_Array();
  a.Init_Array(imin,imax,jmin,jmax,kmin,kmax,halo_size);
  input.read(reinterpret_cast<char *>(a.Raw_Array_Pointer()), sizeof(T)*size);  
} 
//*****************************************************************************
// Helper: Outputs local portion of global VECTOR array to the ofstream
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Read_Binary_Local_Array(ifstream& input, 
    ARRAY_3D<VECTOR_3D<T> >& va)
{
  int size, imin,imax,jmin,jmax,kmin,kmax,halo_size;
  input.read(reinterpret_cast<char *>(&size),sizeof(int));
  input.read(reinterpret_cast<char *>(&imin),sizeof(int));
  input.read(reinterpret_cast<char *>(&imax),sizeof(int));
  input.read(reinterpret_cast<char *>(&jmin),sizeof(int));
  input.read(reinterpret_cast<char *>(&jmax),sizeof(int));
  input.read(reinterpret_cast<char *>(&kmin),sizeof(int));
  input.read(reinterpret_cast<char *>(&kmax),sizeof(int));
  input.read(reinterpret_cast<char *>(&halo_size),sizeof(int));
  va.Delete_Array();
  va.Init_Array(imin,imax,jmin,jmax,kmin,kmax,halo_size);
  T vx[size], vy[size], vz[size];
  input.read(reinterpret_cast<char *>(&vx), sizeof(T)*size); 
  input.read(reinterpret_cast<char *>(&vy), sizeof(T)*size);
  input.read(reinterpret_cast<char *>(&vz), sizeof(T)*size);
  for(int n=0; n<va.Total_Size_With_Halo(); n++) 
    va.Raw_Array(n) = VECTOR_3D<T>(vx[n],vy[n],vz[n]);
} 
//*****************************************************************************
// Writes local T* array to a file '$array_name{_t##}.m' on disk
//*****************************************************************************
template <class T> int MPI_DRIVER<T>::Write_Local_Array_To_Disk(
    string a_name, T* a, int size, int timestep)
{
  // create filename
  std::stringstream filename;
  filename << output_dir;
  filename << a_name;
  if(timestep)
    filename << "_t" << timestep;
  filename << ".m";
  // open file
  ofstream output_file(filename.str().c_str(), ios::out);
  if(!output_file){
    cout<<"ERROR: could not open file for writing"<<endl;
    return 0;
  }
  // write into file
  output_file << a_name;
  if(timestep) output_file <<"_"<<timestep;
  output_file << " = [";
  output_file.precision(15);
  for(int i = 0; i < size; i++) 
    if(i < size-1) output_file<<a[i]<<", "; else output_file<<a[i]<<"];";
  output_file.close(); 
  return 1;
}
/*
//*****************************************************************************
// Exchange ghost values in each dimension for a scalar 3D array.
// If halo=2, then 4 indexes:(-1, 0, n+1, n+2) are exchanged.
// Rules: (n+1,n+2)<-(1,2) and (n-1,n)->(-1,0).
//*****************************************************************************
template <class T>
void MPI_DRIVER<T>::Exchange_Ghost_Values_For_Scalar_Field(ARRAY_3D<T>& scalar,
int Halo_size)
{
if(Halo_size==-1) Halo_size = scalar.Halo_Size(); //default values is '-1'
int min_X=scalar.I_Min() - Halo_size, max_X=scalar.I_Max() + Halo_size,  
min_Y=scalar.J_Min() - Halo_size, max_Y=scalar.J_Max() + Halo_size, 
min_Z=scalar.K_Min() - Halo_size, max_Z=scalar.K_Max() + Halo_size,
X_size=scalar.I_Size()+ 2*Halo_size, Y_size=scalar.J_Size()+ 2*Halo_size,
Z_size=scalar.K_Size()+ 2*Halo_size;

MPI_Request request[4];
MPI_Status status[4];
int request_counter=0;

T recv_west_message[Y_size][Halo_size], send_west_message[Y_size][Halo_size];
T recv_east_message[Y_size][Halo_size], send_east_message[Y_size][Halo_size];
T recv_suth_message[X_size][Halo_size], send_suth_message[X_size][Halo_size];
T recv_nrth_message[X_size][Halo_size], send_nrth_message[X_size][Halo_size];
T recv_back_message[X_size][Halo_size], send_back_message[X_size][Halo_size];
T recv_frnt_message[X_size][Halo_size], send_frnt_message[X_size][Halo_size];

// exchange in I-DIRECTION: east-west

for(int k = min_Z; k <= max_Z; k++) {
//send-receive ghost values  
if(west_proc != MPI_PROC_NULL) {
MPI_Irecv(recv_west_message, Y_size*Halo_size, MPI_DOUBLE, 
west_proc, 0, grid_comm, &request[request_counter++]);
for(int j=min_Y; j <= max_Y; j++)
for(int h=0; h < Halo_size; h++)
send_west_message[j-min_Y][h] = scalar(h+1,j,k); // i=1,2 (halo=2)
MPI_Isend(send_west_message, Y_size*Halo_size, MPI_DOUBLE, 
west_proc, 1, grid_comm, &request[request_counter++]);
}

if(east_proc != MPI_PROC_NULL) {
MPI_Irecv(recv_east_message, Y_size*Halo_size, MPI_DOUBLE, 
east_proc, 1, grid_comm, &request[request_counter++]);
for(int j=min_Y; j <= max_Y; j++)
for(int h=0; h < Halo_size; h++){
int i_index = max_X - 2*Halo_size + h+1; // i=n-1,n (if halo=2)
send_east_message[j-min_Y][h]=scalar(i_index,j,k);
}
MPI_Isend(send_east_message, Y_size*Halo_size, MPI_DOUBLE, 
east_proc, 0, grid_comm, &request[request_counter++]);
}

// wait until exchange completes
MPI_Waitall(request_counter, request, status); 
request_counter = 0;

// assign received values
if(west_proc != MPI_PROC_NULL)
for(int j=min_Y; j <= max_Y; j++)
for(int h=0; h < Halo_size; h++){
int i_index = -Halo_size + (h+1); // i=-1,0 (if halo=2)
scalar(i_index,j,k) = recv_west_message[j-min_Y][h];
}
if(east_proc != MPI_PROC_NULL)
for(int j=min_Y; j <= max_Y; j++)
for(int h=0; h < Halo_size; h++){
int i_index = max_X - Halo_size + (h+1); // i=n+1,n+2 (if halo=2)
scalar(i_index,j,k) = recv_east_message[j-min_Y][h]; 
}

}//for: k

// Check if the exchange succeeded
bool success = true;
for(int k = min_Z; k <= max_Z; k++) 
for(int j= min_Y; j <= max_Y; j++)
if( scalar(-1,j,k) != scalar(max_X-3,j,k) ||
    scalar( 0,j,k) != scalar(max_X-2,j,k) ||
    scalar( 1,j,k) != scalar(max_X-1,j,k) ||
    scalar( 2,j,k) != scalar(max_X,j,k) )
success = false;
if(!success) cout<<"Exchange failed!"<<endl;

// exchange in J-DIRECTION: south-north

for(int k = min_Z; k <= max_Z; k++) {
  //send-receive ghost values  
  if(suth_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_suth_message, X_size*Halo_size, MPI_DOUBLE, 
        suth_proc, 0, grid_comm, &request[request_counter++]);
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++)
        send_suth_message[i-min_X][h] = scalar(i,h+1,k); // j=1,2 (halo=2)
    MPI_Isend(send_suth_message, X_size*Halo_size, MPI_DOUBLE, 
        suth_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(nrth_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_nrth_message, X_size*Halo_size, MPI_DOUBLE, 
        nrth_proc, 1, grid_comm, &request[request_counter++]);
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++){
        int j_index = max_Y - 2*Halo_size + h+1; // j=n-1,n (halo=2)
        send_nrth_message[i-min_X][h]=scalar(i,j_index,k);
      }
    MPI_Isend(send_nrth_message, X_size*Halo_size, MPI_DOUBLE, 
        nrth_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(suth_proc != MPI_PROC_NULL)
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++){
        int j_index = -Halo_size + (h+1); // j=-1,0 (if halo=2)
        scalar(i,j_index,k) = recv_suth_message[i-min_X][h];
      }
  if(nrth_proc != MPI_PROC_NULL)
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++){
        int j_index = max_Y - Halo_size + (h+1); // j=n+1,n+2 (if halo=2)
        scalar(i,j_index,k) = recv_nrth_message[i-min_X][h]; 
      }

}//for: k

// exchange in K-DIRECTION: back-frnt

for(int j = min_Y; j <= max_Y; j++) {

  //send-receive ghost values  
  if(back_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_back_message, X_size*Halo_size, MPI_DOUBLE, 
        back_proc, 0, grid_comm, &request[request_counter++]);
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++)
        send_back_message[i-min_X][h] = scalar(i,j,h+1); // k=1,2 (halo=2)
    MPI_Isend(send_back_message, X_size*Halo_size, MPI_DOUBLE, 
        back_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(frnt_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_frnt_message, X_size*Halo_size, MPI_DOUBLE, 
        frnt_proc, 1, grid_comm, &request[request_counter++]);
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++){
        int k_index = max_Z - 2*Halo_size + h+1; // k=n-1,n (if halo=2)
        send_frnt_message[i-min_X][h]=scalar(i,j,k_index);
      }
    MPI_Isend(send_frnt_message, X_size*Halo_size, MPI_DOUBLE, 
        frnt_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(back_proc != MPI_PROC_NULL)
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++){
        int k_index = -Halo_size+(h+1); // k=-1,0 (if halo=2)
        scalar(i,j,k_index) = recv_back_message[i-min_X][h];
      }
  if(frnt_proc != MPI_PROC_NULL)
    for(int i=min_X; i <= max_X; i++)
      for(int h=0; h < Halo_size; h++){
        int k_index = max_Z - Halo_size + (h+1); // k=n+1,n+2 (if halo=2)
        scalar(i,j,k_index) = recv_frnt_message[i-min_X][h]; 
      }

}//for: j
}
*/
//*****************************************************************************
// Exchange ghost values in each dimension for a scalar 3D array.
// If halo=2, then 4 indexes:(-1, 0, n+1, n+2) are exchanged.
// Rules: (n+1,n+2)<-(1,2) and (n-1,n)->(-1,0).
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Exchange_Ghost_Values_For_Scalar_Field(ARRAY_3D<T>& scalar,
    int Halo_size)
{
  if(Halo_size==-1) Halo_size = scalar.Halo_Size(); //default values is '-1'

  int min_X=scalar.I_Min() - Halo_size, max_X=scalar.I_Max() + Halo_size,  
      min_Y=scalar.J_Min() - Halo_size, max_Y=scalar.J_Max() + Halo_size, 
      min_Z=scalar.K_Min() - Halo_size, max_Z=scalar.K_Max() + Halo_size,
      X_size=scalar.I_Size()+ 2*Halo_size, Y_size=scalar.J_Size()+ 2*Halo_size,
      Z_size=scalar.K_Size()+ 2*Halo_size;

  MPI_Request request[4];
  MPI_Status status[4];
  int request_counter=0;

  T recv_west_message[Z_size][Y_size][Halo_size], 
    send_west_message[Z_size][Y_size][Halo_size],
    recv_east_message[Z_size][Y_size][Halo_size], 
    send_east_message[Z_size][Y_size][Halo_size],
    recv_suth_message[Z_size][X_size][Halo_size], 
    send_suth_message[Z_size][X_size][Halo_size],
    recv_nrth_message[Z_size][X_size][Halo_size], 
    send_nrth_message[Z_size][X_size][Halo_size],
    recv_back_message[Y_size][X_size][Halo_size], 
    send_back_message[Y_size][X_size][Halo_size],
    recv_frnt_message[Y_size][X_size][Halo_size], 
    send_frnt_message[Y_size][X_size][Halo_size];

  // exchange in I-DIRECTION: east-west

  //send-receive ghost values  
  if(west_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_west_message, Z_size*Y_size*Halo_size, MPI_DOUBLE, 
        west_proc, 0, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++)
          send_west_message[k-min_Z][j-min_Y][h] = scalar(h+1,j,k); //i=1,2(h=2)
    MPI_Isend(send_west_message, Z_size*Y_size*Halo_size, MPI_DOUBLE, 
        west_proc, 1, grid_comm, &request[request_counter++]);
  }
  if(east_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_east_message, Z_size*Y_size*Halo_size, MPI_DOUBLE, 
        east_proc, 1, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          int i_index = max_X - 2*Halo_size + h+1; // i=n-1,n (if halo=2)
          send_east_message[k-min_Z][j-min_Y][h]=scalar(i_index,j,k);
        }
    MPI_Isend(send_east_message, Z_size*Y_size*Halo_size, MPI_DOUBLE, 
        east_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(west_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          int i_index = -Halo_size + (h+1); // i=-1,0 (if halo=2)
          scalar(i_index,j,k) = recv_west_message[k-min_Z][j-min_Y][h];
        }
  if(east_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          int i_index = max_X - Halo_size + (h+1); // i=n+1,n+2 (if halo=2)
          scalar(i_index,j,k) = recv_east_message[k-min_Z][j-min_Y][h]; 
        }

  // exchange in J-DIRECTION: south-north

  //send-receive ghost values  
  if(suth_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_suth_message, Z_size*X_size*Halo_size, MPI_DOUBLE, 
        suth_proc, 0, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++)
          send_suth_message[k-min_Z][i-min_X][h] = scalar(i,h+1,k); //j=1,2(h=2)
    MPI_Isend(send_suth_message, Z_size*X_size*Halo_size, MPI_DOUBLE, 
        suth_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(nrth_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_nrth_message, Z_size*X_size*Halo_size, MPI_DOUBLE, 
        nrth_proc, 1, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int j_index = max_Y - 2*Halo_size + h+1; // j=n-1,n (halo=2)
          send_nrth_message[k-min_Z][i-min_X][h]=scalar(i,j_index,k);
        }
    MPI_Isend(send_nrth_message, Z_size*X_size*Halo_size, MPI_DOUBLE, 
        nrth_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(suth_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int j_index = -Halo_size + (h+1); // j=-1,0 (if halo=2)
          scalar(i,j_index,k) = recv_suth_message[k-min_Z][i-min_X][h];
        }
  if(nrth_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int j_index = max_Y - Halo_size + (h+1); // j=n+1,n+2 (if halo=2)
          scalar(i,j_index,k) = recv_nrth_message[k-min_Z][i-min_X][h]; 
        }

  // exchange in K-DIRECTION: back-frnt

  //send-receive ghost values  
  if(back_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_back_message, Y_size*X_size*Halo_size, MPI_DOUBLE, 
        back_proc, 0, grid_comm, &request[request_counter++]);
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++)
          send_back_message[j-min_Y][i-min_X][h] = scalar(i,j,h+1); //k=1,2(h=2)
    MPI_Isend(send_back_message, Y_size*X_size*Halo_size, MPI_DOUBLE, 
        back_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(frnt_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_frnt_message, Y_size*X_size*Halo_size, MPI_DOUBLE, 
        frnt_proc, 1, grid_comm, &request[request_counter++]);
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int k_index = max_Z - 2*Halo_size + h+1; // k=n-1,n (if halo=2)
          send_frnt_message[j-min_Y][i-min_X][h]=scalar(i,j,k_index);
        }
    MPI_Isend(send_frnt_message, Y_size*X_size*Halo_size, MPI_DOUBLE, 
        frnt_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(back_proc != MPI_PROC_NULL)
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int k_index = -Halo_size+(h+1); // k=-1,0 (if halo=2)
          scalar(i,j,k_index) = recv_back_message[j-min_Y][i-min_X][h];
        }
  if(frnt_proc != MPI_PROC_NULL)
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int k_index = max_Z - Halo_size + (h+1); // k=n+1,n+2 (if halo=2)
          scalar(i,j,k_index) = recv_frnt_message[j-min_Y][i-min_X][h]; 
        }
}
//*****************************************************************************
// Exchange ghost values in each dimension for a VECTOR_3D array.
// If halo=2, then 4 indexes:(-1, 0, n+1, n+2) are exchanged.
// Rules: (n+1,n+2)<-(1,2) and (n-1,n)->(-1,0).
//*****************************************************************************
  template <class T>
void MPI_DRIVER<T>::Exchange_Ghost_Values_For_Vector_Field(
    ARRAY_3D<VECTOR_3D<T> >& vector)
{
  int min_X=vector.I_Min_With_Halo(), max_X=vector.I_Max_With_Halo(),  
      min_Y=vector.J_Min_With_Halo(), max_Y=vector.J_Max_With_Halo(), 
      min_Z=vector.K_Min_With_Halo(), max_Z=vector.K_Max_With_Halo(),
      X_size=vector.I_Size_With_Halo(), Y_size=vector.J_Size_With_Halo(),
      Z_size=vector.K_Size_With_Halo(), Halo_size=vector.Halo_Size();
  MPI_Request request[4];
  MPI_Status status[4];
  int request_counter=0;

  T recv_west_message[Z_size][Y_size][Halo_size][3], 
    send_west_message[Z_size][Y_size][Halo_size][3],
    recv_east_message[Z_size][Y_size][Halo_size][3], 
    send_east_message[Z_size][Y_size][Halo_size][3],
    recv_suth_message[Z_size][X_size][Halo_size][3], 
    send_suth_message[Z_size][X_size][Halo_size][3],
    recv_nrth_message[Z_size][X_size][Halo_size][3], 
    send_nrth_message[Z_size][X_size][Halo_size][3],
    recv_back_message[Y_size][X_size][Halo_size][3], 
    send_back_message[Y_size][X_size][Halo_size][3],
    recv_frnt_message[Y_size][X_size][Halo_size][3], 
    send_frnt_message[Y_size][X_size][Halo_size][3];

  // exchange in I-DIRECTION: east-west

  //send-receive ghost values  
  if(west_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_west_message, Z_size*Y_size*Halo_size*3, MPI_DOUBLE, 
        west_proc, 0, grid_comm, &request[request_counter++]);  
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          send_west_message[k-min_Z][j-min_Y][h][0] = vector(h+1,j,k).x;//i=1,2
          send_west_message[k-min_Z][j-min_Y][h][1] = vector(h+1,j,k).y;
          send_west_message[k-min_Z][j-min_Y][h][2] = vector(h+1,j,k).z;
        }
    MPI_Isend(send_west_message, Z_size*Y_size*Halo_size*3, MPI_DOUBLE, 
        west_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(east_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_east_message, Z_size*Y_size*Halo_size*3, MPI_DOUBLE, 
        east_proc, 1, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          int i_index = max_X - 2*Halo_size + h+1; // i=n-1,n (if halo=2)
          send_east_message[k-min_Z][j-min_Y][h][0] = vector(i_index,j,k).x;
          send_east_message[k-min_Z][j-min_Y][h][1] = vector(i_index,j,k).y;
          send_east_message[k-min_Z][j-min_Y][h][2] = vector(i_index,j,k).z;
        }
    MPI_Isend(send_east_message, Z_size*Y_size*Halo_size*3, MPI_DOUBLE, 
        east_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(west_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          int i_index = -Halo_size + (h+1); // i=-1,0 (if halo=2)
          vector(i_index,j,k).x = recv_west_message[k-min_Z][j-min_Y][h][0];
          vector(i_index,j,k).y = recv_west_message[k-min_Z][j-min_Y][h][1];
          vector(i_index,j,k).z = recv_west_message[k-min_Z][j-min_Y][h][2];
        }
  if(east_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int j=min_Y; j <= max_Y; j++)
        for(int h=0; h < Halo_size; h++){
          int i_index = max_X - Halo_size + (h+1); // i=n+1,n+2 (if halo=2)
          vector(i_index,j,k).x = recv_east_message[k-min_Z][j-min_Y][h][0];
          vector(i_index,j,k).y = recv_east_message[k-min_Z][j-min_Y][h][1];
          vector(i_index,j,k).z = recv_east_message[k-min_Z][j-min_Y][h][2]; 
        }

  // exchange in J-DIRECTION: south-north

  //send-receive ghost values  
  if(suth_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_suth_message, Z_size*X_size*Halo_size*3, MPI_DOUBLE, 
        suth_proc, 0, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          send_suth_message[k-min_Z][i-min_X][h][0] = vector(i,h+1,k).x; //j=1,2
          send_suth_message[k-min_Z][i-min_X][h][1] = vector(i,h+1,k).y;
          send_suth_message[k-min_Z][i-min_X][h][2] = vector(i,h+1,k).z;
        }
    MPI_Isend(send_suth_message, Z_size*X_size*Halo_size*3, MPI_DOUBLE, 
        suth_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(nrth_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_nrth_message, Z_size*X_size*Halo_size*3, MPI_DOUBLE, 
        nrth_proc, 1, grid_comm, &request[request_counter++]);
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int j_index = max_Y - 2*Halo_size + h+1; // j=n-1,n (halo=2)
          send_nrth_message[k-min_Z][i-min_X][h][0] = vector(i,j_index,k).x;
          send_nrth_message[k-min_Z][i-min_X][h][1] = vector(i,j_index,k).y;
          send_nrth_message[k-min_Z][i-min_X][h][2] = vector(i,j_index,k).z;
        }
    MPI_Isend(send_nrth_message, Z_size*X_size*Halo_size*3, MPI_DOUBLE, 
        nrth_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(suth_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int j_index = -Halo_size + (h+1); // j=-1,0 (if halo=2)
          vector(i,j_index,k).x = recv_suth_message[k-min_Z][i-min_X][h][0];
          vector(i,j_index,k).y = recv_suth_message[k-min_Z][i-min_X][h][1];
          vector(i,j_index,k).z = recv_suth_message[k-min_Z][i-min_X][h][2];
        }
  if(nrth_proc != MPI_PROC_NULL)
    for(int k=min_Z; k <= max_Z; k++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int j_index = max_Y - Halo_size + (h+1); // j=n+1,n+2 (if halo=2)
          vector(i,j_index,k).x = recv_nrth_message[k-min_Z][i-min_X][h][0]; 
          vector(i,j_index,k).y = recv_nrth_message[k-min_Z][i-min_X][h][1]; 
          vector(i,j_index,k).z = recv_nrth_message[k-min_Z][i-min_X][h][2]; 
        }

  // exchange in K-DIRECTION: back-frnt

  //send-receive ghost values  
  if(back_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_back_message, Y_size*X_size*Halo_size*3, MPI_DOUBLE, 
        back_proc, 0, grid_comm, &request[request_counter++]);
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          send_back_message[j-min_Y][i-min_X][h][0] = vector(i,j,h+1).x;//k=1,2
          send_back_message[j-min_Y][i-min_X][h][1] = vector(i,j,h+1).y;
          send_back_message[j-min_Y][i-min_X][h][2] = vector(i,j,h+1).z;
        }
    MPI_Isend(send_back_message, Y_size*X_size*Halo_size*3, MPI_DOUBLE, 
        back_proc, 1, grid_comm, &request[request_counter++]);
  }

  if(frnt_proc != MPI_PROC_NULL) {
    MPI_Irecv(recv_frnt_message, Y_size*X_size*Halo_size*3, MPI_DOUBLE, 
        frnt_proc, 1, grid_comm, &request[request_counter++]);
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int k_index = max_Z - 2*Halo_size + h+1; // k=n-1,n (if halo=2)
          send_frnt_message[j-min_Y][i-min_X][h][0] = vector(i,j,k_index).x;
          send_frnt_message[j-min_Y][i-min_X][h][1] = vector(i,j,k_index).y;
          send_frnt_message[j-min_Y][i-min_X][h][2] = vector(i,j,k_index).z;
        }
    MPI_Isend(send_frnt_message, Y_size*X_size*Halo_size*3, MPI_DOUBLE, 
        frnt_proc, 0, grid_comm, &request[request_counter++]);
  }

  // wait until exchange completes
  MPI_Waitall(request_counter, request, status); 
  request_counter = 0;

  // assign received values
  if(back_proc != MPI_PROC_NULL)
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int k_index = -Halo_size+(h+1); // k=-1,0 (if halo=2)
          vector(i,j,k_index).x = recv_back_message[j-min_Y][i-min_X][h][0];
          vector(i,j,k_index).y = recv_back_message[j-min_Y][i-min_X][h][1];
          vector(i,j,k_index).z = recv_back_message[j-min_Y][i-min_X][h][2];
        }
  if(frnt_proc != MPI_PROC_NULL)
    for(int j=min_Y; j <= max_Y; j++)
      for(int i=min_X; i <= max_X; i++)
        for(int h=0; h < Halo_size; h++){
          int k_index = max_Z - Halo_size + (h+1); // k=n+1,n+2 (if halo=2)
          vector(i,j,k_index).x = recv_frnt_message[j-min_Y][i-min_X][h][0];
          vector(i,j,k_index).y = recv_frnt_message[j-min_Y][i-min_X][h][1];
          vector(i,j,k_index).z = recv_frnt_message[j-min_Y][i-min_X][h][2];
        }
}
//*****************************************************************************
// Finds MAX value among all procs and replaces 'value' with MAX for everyone.
// Note: used in CFL calculations.
//*****************************************************************************
template<class T>
void MPI_DRIVER<T>::Replace_With_Max_Value_Among_All_Procs(T& value){
  T max_value = (T)0;
  //find max value (on proc==0) among all CPUs and broadcast it to everyone
  MPI_Reduce(&value, &max_value, 1, MPI_DOUBLE, MPI_MAX, 0, grid_comm);
  MPI_Bcast(&max_value, 1, MPI_DOUBLE, 0, grid_comm);
  value = max_value;
}
//*****************************************************************************
// Integer version of previous function
//*****************************************************************************
template<class T>
void MPI_DRIVER<T>::Replace_With_Max_Value_Among_All_Procs(int& value){
  int max_value = 0;
  //find max value (on proc==0) among all CPUs and broadcast it to everyone
  MPI_Reduce(&value, &max_value, 1, MPI_INT, MPI_MAX, 0, grid_comm);
  MPI_Bcast(&max_value, 1, MPI_INT, 0, grid_comm);
  value = max_value;
}
//*****************************************************************************
// Finds SUM of values and brodcasts it to all procs.
// Note: used in potential background energy compuation.
//*****************************************************************************
template<class T>
void MPI_DRIVER<T>::Replace_With_Sum_On_All_Procs(T& value){
  T sum = (T)0;
  //collect sum of values on proc==0 broadcast it to everyone
  MPI_Reduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm);
  MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, grid_comm);
  value = sum;
}
//*****************************************************************************
// Helper function to check if the output directory is missing and create it
//*****************************************************************************
template<class T>
void MPI_DRIVER<T>::Create_Output_Directory_If_Missing(){
  //check if output directory exists
  struct stat st;
  if(stat(output_dir.c_str(),&st) != 0){
    cout<<"ERROR: output directory does not exist"<<endl;
    cout<<"Creating './output' directory...";
    int status = mkdir(output_dir.c_str(), S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    if(status==0) cout<<"SUCCESS"<<endl; else cout<<"FAILURE"<<endl;
  }
}
//*****************************************************************************
#endif
