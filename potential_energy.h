// Potential energy class
#ifndef __POTENTIAL_ENERGY__
#define __POTENTIAL_ENERGY__

#include "curvilinear_grid.h"

template<class T=double>
struct CELL
{
  T rho, volume, z, z_star, local_height, laplacian_rho;
};

template<class T=double>
class POTENTIAL_ENERGY
{
 public:
  POTENTIAL_ENERGY(PARAMETERS<T> *pa, MPI_DRIVER<T> *md, CURVILINEAR_GRID<T> *g,
                   ARRAY_3D<T> *rhoin) 
 : parameters(pa), mpi_driver(md), grid(g), rho(rhoin), p(md->total_procs)
  {
    E_background = new T[parameters->max_timestep]; 
    E_potential = new T[parameters->max_timestep]; 
    Phi_d = new T[parameters->max_timestep];
    F_background = new T[parameters->max_timestep]; 
    F_potential = new T[parameters->max_timestep]; 
    array_size = 0;
    local_sorted_array_size = 0;
    rho_sorted_cells = NULL;
    local_array_size = rho->Total_Size();
    //find size of global 'rho' array
    MPI_Reduce(&local_array_size, &global_array_size, 1, MPI_INT, 
	       MPI_SUM, 0, md->Get_Grid_Communicator());
    MPI_Bcast(&global_array_size, 1, MPI_INT, 0,md->Get_Grid_Communicator());
    rho_local_cells = new CELL<T>[local_array_size];

    rho_local_samples = new CELL<T>[p];
    if(!mpi_driver->my_rank) 
      rho_global_samples = new CELL<T>[p*p]; else rho_global_samples=NULL;
    pivots = new CELL<T>[p-1];
    sorted_array_part_sizes_send = new int[p]; 
    sorted_array_part_sizes_recv = new int[p];  
    send_displacements = new int[p]; 
    recv_displacements = new int[p];
    // register MPI datatype
    int block_count[1] = {6}; //old 2
    MPI_Aint offset[1] = {0};
    MPI_Datatype old_types[1] = {MPI_DOUBLE};
    MPI_Type_struct(1, block_count, offset, old_types, &cell_type);
    MPI_Type_commit(&cell_type);
  }

  ~POTENTIAL_ENERGY() 
  {
    delete[] E_background; delete[] E_potential;
    delete[] F_background; delete[] F_potential; delete[] Phi_d;
    delete[] rho_local_cells; delete[] rho_local_samples;
    if(!mpi_driver->my_rank) delete[] rho_global_samples;
    delete[] pivots;
    delete[] sorted_array_part_sizes_send;delete[] sorted_array_part_sizes_recv;
    delete[] send_displacements; delete[] recv_displacements;
    MPI_Type_free(&cell_type);
  }

  T Calculate(){
    E_background[array_size] = Background_Potential_Energy();
    E_potential[array_size] = Total_Potential_Energy();
    Phi_d[array_size] = Mixing();
    if(parameters->west_velocity){
      F_background[array_size] = Background_Potential_Energy_Flux();
      F_potential[array_size] = Potential_Energy_Flux();
    }
    return E_background[array_size++];
  }

  int Write_To_Disk();

 private:
  T Background_Potential_Energy();
  T Background_Potential_Energy_Flux();
  T Mixing();
  T Total_Potential_Energy();
  T Potential_Energy_Flux();

  void Convert_ARRAY_3D_To_Linear_Array();
  void Sort_Global_Density_Array();
  void Aggregate_All_Samples_On_Root();
  void Pick_Pivots_From_Global_Samples();
  void Split_Arrays_Based_On_Pivots();
  void Redistribute_Local_Arrays();

  T Calculate_Planform_Area(T cell_height);
  T Calculate_Planform_Area_Smooth_Grid(T cell_height, T slope, T xend, T zend, 
    T radius, T xc, T zc, T x_length, T y_length, T z_length);
  T Calculate_Cell_Height(T cell_volume, T cell_bottom_height,
    T slope, T x_slope, T x_length, T y_length);
  T Calculate_Cell_Centroid_Height(T local_height, T cell_bottom_height,
    T slope, T x_slope, T x_length, T y_length);
  T Calculate_Laplacian_of_Rho(int i, int j, int k);
  T Receive_Initial_Local_Height();
  void Send_Final_Local_Height(T final_cell_height);

  static int compare_cells(const void *a, const void *b);

  CURVILINEAR_GRID<T> *grid;
  PARAMETERS<T> *parameters;
  MPI_DRIVER<T> *mpi_driver;
  ARRAY_3D<T> *rho;
  int &p; //shortcut for num cpu
  int local_array_size, global_array_size, local_sorted_array_size;
  CELL<T> *rho_local_cells, *rho_sorted_cells;
  CELL<T> *rho_local_samples, *rho_global_samples, *pivots;
  int *sorted_array_part_sizes_send, *sorted_array_part_sizes_recv,
      *send_displacements, *recv_displacements;
  int array_size;
  T *E_background, *E_potential, *Phi_d, *F_background, *F_potential;
  T *sorted_z_star, *sorted_local_height, *sorted_laplacian_rho;
  MPI_Datatype cell_type;
};
//*****************************************************************************
// Background Potential Energy: measure of numerical diffusion
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Background_Potential_Energy()
{
  /*
  //For rectangular domain

  T E_b = (T)0;
  Convert_ARRAY_3D_To_Linear_Array(); 
  Sort_Global_Density_Array();
 
  T cell_height = Receive_Initial_Local_Height(); 
  inv_domain_planform_width = (T)1 / (parameters->x_length*parameters->y_length); 
  //cout<<"Initial Height="<<cell_height<<" on CPU#"<<mpi_driver->my_rank<<endl;
  //cout<<"Sorted array size="<<local_sorted_array_size<<" out of "
  //    <<local_array_size <<" on CPU#"<<mpi_driver->my_rank<<endl;

  // calculate local E_b
  for(int n = 0; n < local_sorted_array_size; n++){
    T local_height = rho_sorted_cells[n].volume * inv_domain_planform_width;
    if(!mpi_driver->my_rank && !n) local_height *= (T).5; //first cell
    cell_height += local_height;
    rho_sorted_cells[n].z_star = cell_height;
    E_b += rho_sorted_cells[n].rho * rho_sorted_cells[n].volume * cell_height;
    //if(n && rho_sorted_cells[n].rho > rho_sorted_cells[n-1].rho){
    //  cout.precision(15);
    //  cout<<"NOT SORTED:"<<n<<":"<< rho_sorted_cells[n].rho 
    //      <<" > "<<rho_sorted_cells[n-1].rho<<endl;
    //}
  }
  E_b *= parameters->g;
  Send_Final_Local_Height(cell_height);
  //cout<<"E_b="<<E_b<<" on CPU#"<<mpi_driver->my_rank<<endl;
  // sum over all procs
  if(p>1) {
    mpi_driver->Replace_With_Sum_On_All_Procs(E_b);
  }
  return E_b;
  */

  //For sloping domain 
  T E_b = (T)0;
  T local_height, centroid_height;
  Convert_ARRAY_3D_To_Linear_Array(); 
  Sort_Global_Density_Array();
  T cell_height = Receive_Initial_Local_Height(); 
  
  // calculate local E_b
  for(int n = 0; n < local_sorted_array_size; n++){
    local_height = Calculate_Cell_Height(rho_sorted_cells[n].volume, 
      cell_height,parameters->slope,parameters->x_s,parameters->x_length,parameters->y_length);
    rho_sorted_cells[n].local_height = local_height;
    centroid_height = Calculate_Cell_Centroid_Height(local_height,cell_height,
      parameters->slope,parameters->x_s,parameters->x_length,parameters->y_length);
    rho_sorted_cells[n].z_star = centroid_height;
    E_b += rho_sorted_cells[n].rho * rho_sorted_cells[n].volume 
         * centroid_height;
    cell_height += local_height;
    }
  E_b *= parameters->g;
  Send_Final_Local_Height(cell_height);

  // sum over all procs
  if(p>1) mpi_driver->Replace_With_Sum_On_All_Procs(E_b);

  return E_b;
}
//*****************************************************************************
// Calculate phi_d, irreversible diapycnal mixing
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Mixing()
{
  T phi_d = (T)0;

  // calculate Phi_d
  for(int n = 0; n < local_sorted_array_size; n++){
    phi_d += rho_sorted_cells[n].z_star * rho_sorted_cells[n].laplacian_rho 
           * rho_sorted_cells[n].volume;
  }
  phi_d *= parameters->g * parameters->molecular_diffusivity;

  // sum over all procs
  if(p>1) {
    mpi_driver->Replace_With_Sum_On_All_Procs(phi_d);
    if(!parameters->west_velocity)
      if(rho_sorted_cells){ 
        delete[] rho_sorted_cells; //created in Sorting func
      }
  }

  return phi_d;
}
//*****************************************************************************
// Background Potential Energy Flux at West Boundary
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Background_Potential_Energy_Flux()
{
  T F_Eb = (T)0;
  
  // calculate F_Eb
  for(int n = 0; n < local_sorted_array_size; n++){
    F_Eb += rho_sorted_cells[n].rho  
          * rho_sorted_cells[n].z_star 
          * (rho_sorted_cells[n].local_height*parameters->y_length)
          * parameters->forcing_amp*cos(parameters->m
                                       *(rho_sorted_cells[n].z_star))
                                   *sin(parameters->freq*parameters->time);
  }
  F_Eb *= parameters->g;

  // sum over all procs
  if(p>1) {
    mpi_driver->Replace_With_Sum_On_All_Procs(F_Eb);
    if(rho_sorted_cells){ 
      delete[] rho_sorted_cells; //created in Sorting func
    }
  }

  return F_Eb;
}
//*****************************************************************************
// Calculate planform area of domain for E_b calculation
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Calculate_Planform_Area(T cell_height)
{
  T area = (T)0;
    
  if(cell_height < (parameters->slope * (parameters->x_length - parameters->x_s)))
    area = (parameters->x_s + (1/parameters->slope)*cell_height);
  else
    area = parameters->x_length;

  area *= parameters->y_length;
  return area;
}
//*****************************************************************************
// Calculate planform area of smooth grid domain for E_b calculation
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Calculate_Planform_Area_Smooth_Grid(T cell_height, T slope, T xend, T zend, 
    T radius, T xc, T zc, T x_length, T y_length, T z_length)
{
  T area = (T)0;
    
  if(cell_height < slope*x_length + (zend - slope*xend) + z_length && cell_height > zend + z_length)
    area = ((cell_height - z_length) - (zend - slope*xend))/slope; 
  else if(cell_height <= zend + z_length)
    area = sqrt(pow(radius,2) - pow((cell_height - z_length)-zc,2)) + xc; 
  else
    area = parameters->x_length;

  area *= parameters->y_length;
  return area;
}
//*****************************************************************************
// Calculate local cell height for Eb calculation
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Calculate_Cell_Height(T cell_volume, T cell_bottom_height,
    T slope, T x_slope, T x_length, T y_length)
{
  T local_height, A_bottom, extra_volume, extra_height;

  if(cell_bottom_height < (slope * (x_length - x_slope))){
    A_bottom = y_length * (x_slope + (1/slope)*cell_bottom_height);
    local_height = (slope/y_length) 
                 * (-A_bottom + sqrt(pow(A_bottom,2) + 2*cell_volume*y_length/slope));
    if((cell_bottom_height + local_height) > (slope * (x_length - x_slope))){
      extra_height = (cell_bottom_height + local_height) 
                   - (slope * (x_length - x_slope));
      extra_volume = 0.5*y_length*pow(extra_height,2)/slope;
      local_height += extra_volume / (y_length * x_length);
    }
  }
  else
    local_height = cell_volume / (y_length * x_length); 

  return local_height;
}
//*****************************************************************************
// Calculate cell centroid height for Eb calculation
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Calculate_Cell_Centroid_Height(T local_height, T cell_bottom_height,
    T slope, T x_slope, T x_length, T y_length)
{
  T centroid_height, a, b;

  if(cell_bottom_height < (slope * (x_length - x_slope))){
    a = x_slope + (1/slope)*cell_bottom_height;
    b = x_slope + (1/slope)*(cell_bottom_height + local_height);
    centroid_height = (local_height/3)*((2*a+b)/(a+b));
  }
  else
    centroid_height = 0.5*local_height;

  centroid_height += cell_bottom_height;
  return centroid_height;
}
//*****************************************************************************
// Total Potential Energy
//*****************************************************************************
  template<class T> 
T POTENTIAL_ENERGY<T>::Total_Potential_Energy()
{
  T E_p = (T)0;
  // calculate local E_p
  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
        T cell_volume = (T)1 / (*grid->inverse_Jacobian)(i,j,k),
          cell_height = (*grid)(i-1,j-1,k-1).z + (*grid)(i+1,j-1,k-1).z
                      + (*grid)(i-1,j+1,k-1).z + (*grid)(i+1,j+1,k-1).z
                      + (*grid)(i-1,j-1,k+1).z + (*grid)(i+1,j-1,k+1).z
                      + (*grid)(i-1,j+1,k+1).z + (*grid)(i+1,j+1,k+1).z;
        cell_height /= (T)8;
        cell_height -= parameters->z_min;
        E_p += (*rho)(i,j,k) * cell_volume * cell_height;
      }
  E_p *= parameters->g;
  // sum over all procs
  if(p>1) mpi_driver->Replace_With_Sum_On_All_Procs(E_p);
  return E_p;
}
//*****************************************************************************
// Potential Energy Flux at West Boundary
//*****************************************************************************
  template<class T> 
T POTENTIAL_ENERGY<T>::Potential_Energy_Flux()
{
  T F_Ep = (T)0;
  // calculate local F_Ep
  if(mpi_driver->west_proc == MPI_PROC_NULL)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
        int i = grid->I_Min();
          T cell_area = ((T)1 / (*grid->inverse_Jacobian)(i,j,k)) 
                      / (parameters->x_length/parameters->num_total_nodes_x),
            cell_height = (*grid)(i-1,j-1,k-1).z + (*grid)(i+1,j-1,k-1).z
                      + (*grid)(i-1,j+1,k-1).z + (*grid)(i+1,j+1,k-1).z
                      + (*grid)(i-1,j-1,k+1).z + (*grid)(i+1,j-1,k+1).z
                      + (*grid)(i-1,j+1,k+1).z + (*grid)(i+1,j+1,k+1).z;
          cell_height /= (T)8;
          cell_height -= parameters->z_min;
          F_Ep += (*rho)(i,j,k) * cell_height * cell_area 
                * (*parameters->west_velocity)(j,k).x;
      }
  F_Ep *= parameters->g;
  // sum over all procs
  if(p>1) mpi_driver->Replace_With_Sum_On_All_Procs(F_Ep);
  return F_Ep;
}
//*****************************************************************************
// convert ARRAY_3D to CELL*
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Convert_ARRAY_3D_To_Linear_Array()
{
  int index = 0;
  for(int i = grid->I_Min(); i <= grid->I_Max(); i++)
    for(int j = grid->J_Min(); j <= grid->J_Max(); j++)
      for(int k = grid->K_Min(); k <= grid->K_Max(); k++){
	      rho_local_cells[index].rho = (*rho)(i,j,k);
        rho_local_cells[index].volume =(T)1/(*grid->inverse_Jacobian)(i,j,k);
        rho_local_cells[index].z = (*grid)(i,j,k).z;
        rho_local_cells[index].z_star = (T)0; //calculated with Eb
        rho_local_cells[index].local_height = (T)0; //calculated with Eb
        rho_local_cells[index++].laplacian_rho = Calculate_Laplacian_of_Rho(i,j,k);
  }
}
//*****************************************************************************
// Sort distributed global density array.
// Return sorted density array with their corresponding cell volumes.
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Calculate_Laplacian_of_Rho(int i, int j, int k)
{
  T laplacian = (T)0;

  //12 and 13
  laplacian +=
      (*grid->G12)(i  ,j,k) * ((*rho)(i  ,j+1,k) - (*rho)(i  ,j-1,k)
                            +  (*rho)(i+1,j+1,k) - (*rho)(i+1,j-1,k))
    - (*grid->G12)(i-1,j,k) * ((*rho)(i  ,j+1,k) - (*rho)(i  ,j-1,k)
                            +  (*rho)(i-1,j+1,k) - (*rho)(i-1,j-1,k));

  laplacian +=
      (*grid->G13)(i  ,j,k) * ((*rho)(i  ,j,k+1) - (*rho)(i  ,j,k-1)
                            +  (*rho)(i+1,j,k+1) - (*rho)(i+1,j,k-1))
    - (*grid->G13)(i-1,j,k) * ((*rho)(i  ,j,k+1) - (*rho)(i  ,j,k+1)
                            +  (*rho)(i-1,j,k+1) - (*rho)(i-1,j,k-1));

  //23 and 21
  laplacian +=
      (*grid->G23)(i,j  ,k) * ((*rho)(i,j  ,k+1) - (*rho)(i,j  ,k-1)
                            +  (*rho)(i,j+1,k+1) - (*rho)(i,j+1,k-1))
    - (*grid->G23)(i,j-1,k) * ((*rho)(i,j  ,k+1) - (*rho)(i,j  ,k-1)
                            +  (*rho)(i,j-1,k+1) - (*rho)(i,j-1,k-1));
  laplacian +=
      (*grid->G21)(i,j  ,k) * ((*rho)(i+1,j  ,k) - (*rho)(i-1,j  ,k)
                            +  (*rho)(i+1,j+1,k) - (*rho)(i-1,j+1,k))
    - (*grid->G21)(i,j-1,k) * ((*rho)(i+1,j  ,k) - (*rho)(i-1,j  ,k)
                            +  (*rho)(i+1,j-1,k) - (*rho)(i-1,j-1,k));

  //31 and 32
  laplacian +=
      (*grid->G31)(i,j,k  ) * ((*rho)(i+1,j,k  ) - (*rho)(i-1,j,k  )
                            +  (*rho)(i+1,j,k+1) - (*rho)(i-1,j,k+1))
    - (*grid->G31)(i,j,k-1) * ((*rho)(i+1,j,k  ) - (*rho)(i-1,j,k  )
                            +  (*rho)(i+1,j,k-1) - (*rho)(i-1,j,k-1));
  laplacian +=
      (*grid->G32)(i,j,k  ) * ((*rho)(i,j+1,k  ) - (*rho)(i,j-1,k  )
                            +  (*rho)(i,j+1,k+1) - (*rho)(i,j-1,k+1))
    - (*grid->G32)(i,j,k-1) * ((*rho)(i,j+1,k  ) - (*rho)(i,j-1,k  )
                            +  (*rho)(i,j+1,k-1) - (*rho)(i,j-1,k-1));
  //11
  laplacian += 
      (*grid->G11)(i  ,j,k) * ((*rho)(i+1,j,k)-(*rho)(i  ,j,k))
    - (*grid->G11)(i-1,j,k) * ((*rho)(i  ,j,k)-(*rho)(i-1,j,k));
  //22
  laplacian +=
      (*grid->G22)(i,j  ,k) * ((*rho)(i,j+1,k)-(*rho)(i,j  ,k))
    - (*grid->G22)(i,j-1,k) * ((*rho)(i,j  ,k)-(*rho)(i,j-1,k));
  //33
  laplacian +=
      (*grid->G33)(i,j,k  ) * ((*rho)(i,j,k+1)-(*rho)(i,j,k  ))
    - (*grid->G33)(i,j,k-1) * ((*rho)(i,j,k  )-(*rho)(i,j,k-1));
  
  laplacian *= (*grid->inverse_Jacobian)(i,j,k);

  return laplacian;
}
//*****************************************************************************
// Sort distributed global density array.
// Return sorted density array with their corresponding cell volumes.
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Sort_Global_Density_Array()
{
  // qsort the local array of cells (rho+volume)
  qsort(rho_local_cells, local_array_size, 
  	                 sizeof(CELL<T>), POTENTIAL_ENERGY<T>::compare_cells);
  if(p>1){
    Aggregate_All_Samples_On_Root();
    Pick_Pivots_From_Global_Samples();
    Split_Arrays_Based_On_Pivots();
    Redistribute_Local_Arrays();
    qsort(rho_sorted_cells, local_sorted_array_size, 
  	                 sizeof(CELL<T>), POTENTIAL_ENERGY<T>::compare_cells);
  }else{
    rho_sorted_cells = rho_local_cells; 
    local_sorted_array_size = local_array_size;
  }
}
//*****************************************************************************
// pick density samples on all Nodes and collect them on Node 0
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Aggregate_All_Samples_On_Root()
{
  for(int i=0; i<p; i++) {
     long int index = (long int) i*(global_array_size/(p*p)); //prevent overflow
     rho_local_samples[i] = rho_local_cells[index];
   }
  MPI_Gather(rho_local_samples,  p, cell_type, 
             rho_global_samples, p, cell_type, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}
//*****************************************************************************
// select pivots out of global samples on Node#0 and brodcast them
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Pick_Pivots_From_Global_Samples()
{  
  if(!mpi_driver->my_rank){
     // sort global_samples
     qsort(rho_global_samples, p*p, sizeof(CELL<T>), 
                                            POTENTIAL_ENERGY<T>::compare_cells);
     // select pivots on Node#0
     for(int i = 0; i < p-1; i++){       
       pivots[i] = rho_global_samples[(i+1)*p + (int)floor(p/2) - 1];
       //cout.precision(15);
       //cout<<"Pivot #"<<i<<"="<<pivots[i].rho<<endl;
     }
  }
  // broadcast pivots from Node#0
  MPI_Bcast(pivots, p-1, cell_type, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}
//*****************************************************************************
// 1) split local_array in p-parts using pivot-separators
// 2) Alltoall exchange of array_part sizes from and to each Node
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Split_Arrays_Based_On_Pivots()
{
  // 1)
  for(int i=0; i<p; i++) sorted_array_part_sizes_send[i] = 0; //reset
  int n, i;
  for(n=0, i=0; n<local_array_size && i<p-1; n++)
    if(rho_local_cells[n].rho >= pivots[i].rho) 
       sorted_array_part_sizes_send[i]++; 
     else 
       {i++;n--;}
   // variable n holds an index +1 after last pivot
   sorted_array_part_sizes_send[p-1] = local_array_size - n;   

   // 2)
   MPI_Alltoall(sorted_array_part_sizes_send, 1, MPI_INT, 
                sorted_array_part_sizes_recv, 1, MPI_INT, MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
}
//*****************************************************************************
// 1) calculate size and allocate space for a new local array
// 2) Alltoall exchange of sub-array parts with known sizes/displacements
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Redistribute_Local_Arrays()
{
  // 1)
  local_sorted_array_size = 0;
  for(int i=0;i<p;i++) local_sorted_array_size+=sorted_array_part_sizes_recv[i];
  rho_sorted_cells = new CELL<T>[local_sorted_array_size];

  // 2)
  send_displacements[0] = 0; recv_displacements[0] = 0;
  for(int i = 1; i < p; i++){ 
    send_displacements[i] = 
                  send_displacements[i-1] + sorted_array_part_sizes_send[i-1];
    recv_displacements[i] = 
                  recv_displacements[i-1] + sorted_array_part_sizes_recv[i-1];
  }
  MPI_Alltoallv(
    rho_local_cells,sorted_array_part_sizes_send,send_displacements, cell_type,
    rho_sorted_cells,sorted_array_part_sizes_recv,recv_displacements, cell_type,
              MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}
//*****************************************************************************
// Receive initial height from the previous proc. In case of #0, domain bottom.
//*****************************************************************************
template<class T> 
T POTENTIAL_ENERGY<T>::Receive_Initial_Local_Height()
{ 
  if(!mpi_driver->my_rank) return (T)0.; //parameters->z_min;
  else{
    T height = (T)0; 
    MPI_Status status; 
    MPI_Recv(&height, 1, MPI_DOUBLE, mpi_driver->my_rank-1, 0, 
	                                               MPI_COMM_WORLD, &status);
    return height;
  }
}
//*****************************************************************************
// Send the top z-value of the current proc (for Receive_Initial_Local_Height())
//*****************************************************************************
template<class T> 
void POTENTIAL_ENERGY<T>::Send_Final_Local_Height(T final_cell_height)
{
  if(mpi_driver->my_rank < mpi_driver->total_procs-1)
    MPI_Send(&final_cell_height, 1, MPI_DOUBLE, mpi_driver->my_rank+1, 0, 
	                                                       MPI_COMM_WORLD); 
}
//*****************************************************************************
// Write to disk up to the current time step from Root
//*****************************************************************************
template<class T> 
int POTENTIAL_ENERGY<T>::Write_To_Disk()
{
  if(!mpi_driver->my_rank){

    stringstream filename;

    //overwrite existing file
    if(parameters->time_step==parameters->save_data_timestep_period){
      filename << parameters->output_dir << "potential_energy.0";
      ofstream output(filename.str().c_str(), ios::out | ios::trunc | ios::binary);
      if(!output){
        cout << "ERROR: could not open potential_energy file for writing" <<endl;
        return 0;
      }
      output.write(reinterpret_cast<char *>(&E_background[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&E_potential[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&Phi_d[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&F_background[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&F_potential[array_size-1]),sizeof(T));
      output.close();
    } 

    //append at end of existing file
    else{
      filename << parameters->output_dir << "potential_energy.0";
      ofstream output(filename.str().c_str(), ios::out | ios::app | ios::binary);
      if(!output){
        cout << "ERROR: could not open potential_energy file for writing" <<endl;
        return 0;
      }
      output.write(reinterpret_cast<char *>(&E_background[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&E_potential[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&Phi_d[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&F_background[array_size-1]),sizeof(T));
      output.write(reinterpret_cast<char *>(&F_potential[array_size-1]),sizeof(T));
      output.close();
    }
    /*
    mpi_driver->Write_Local_Array_To_Disk("E_background", E_background, 
					  array_size, parameters->time_step);
    mpi_driver->Write_Local_Array_To_Disk("E_potential", E_potential, 
					  array_size, parameters->time_step);
    */
  }

  return 1;
}
//*****************************************************************************
// Auxiliary function for qsort(): decreasing order
//*****************************************************************************
template<class T> 
int POTENTIAL_ENERGY<T>::compare_cells(const void *a, const void *b)
{
  if((*(CELL<T>*)a).rho > (*(CELL<T>*)b).rho)
    return -1;
  return (*(CELL<T>*)a).rho < (*(CELL<T>*)b).rho;
}
//*****************************************************************************
#endif
