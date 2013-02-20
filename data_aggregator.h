// Aggregation of time series for physical parameters (e.g., u, P)
#ifndef __DATA_AGGREGATOR__
#define __DATA_AGGREGATOR__

#include "curvilinear_grid.h"

template<class T=double>
class DATA_AGGREGATOR
{
 public:
  DATA_AGGREGATOR(PARAMETERS<T> *pa, MPI_DRIVER<T> *md, CURVILINEAR_GRID<T> *g,
                   ARRAY_3D<VECTOR_3D<T> > *u_in) 
 : parameters(pa), mpi_driver(md), grid(g), u(u_in)
  {

    int host_cpu_on_host = -1;
    array_size = 0;
    // find the host CPU
    bool found = false;
    VECTOR_3D<T> target(.51*parameters->x_length, -.51*parameters->y_length, 
			.51*parameters->z_length);
    for(int i=grid->I_Min(); i<=grid->I_Max()-1 && !found; i++)
      for(int j=grid->J_Min(); j<=grid->J_Max()-1 && !found; j++)
	for(int k=grid->K_Min(); k<=grid->K_Max()-1 && !found; k++){
	  T diam = ((*grid)(i,j,k) - (*grid)(i+1,j+1,k+1)).Magnitude();
	  if( (((*grid)(i,j,k) - target).Magnitude() < diam) &&
              (((*grid)(i+1,j+1,k+1) - target).Magnitude() < diam)){
	    found = true;
	    host_cpu_on_host = mpi_driver->my_rank;
	    center = VECTOR_3D<int>(i,j,k);
	    cout<<"Center on "<<host_cpu_on_host<<" (i,j,k)="<<center
                <<"="<<(*grid)(i,j,k)<<endl;
	  }
	}
    MPI_Reduce(&host_cpu_on_host, &host_cpu, 1, MPI_INT, MPI_MAX, 0, 
	                               md->Get_Grid_Communicator());
    MPI_Bcast(&host_cpu, 1, MPI_INT, 0,md->Get_Grid_Communicator());
    cout<<"On "<<mpi_driver->my_rank<<" host is "<<host_cpu<<endl;
    if(mpi_driver->my_rank==host_cpu)
      u_center_aggregate = new T[parameters->max_timestep];  
    // time averaged velocity array fraction on each CPU
    u_time_average = new ARRAY_3D<VECTOR_3D<T> >(*u);
    u_time_average->Set_All_Elements_To(VECTOR_3D<T>(0.,0.,0.));
    ts_counter = 0;
  }

  ~DATA_AGGREGATOR() {
    if(mpi_driver->my_rank==host_cpu) delete[] u_center_aggregate;
    delete u_time_average;
  }

  void Aggregate(){
    // insert an element
    if(mpi_driver->my_rank==host_cpu) 
      u_center_aggregate[++array_size] = 
	(*u)(center.x,center.y,center.z).Magnitude();
    // add intantaneous V field
    *u_time_average += *u; 
    ts_counter++;
  }

  void Write_To_Disk(){
    if(mpi_driver->my_rank==host_cpu)
    mpi_driver->Write_Local_Array_To_Disk("u_center_aggregate", 
		        u_center_aggregate, array_size, parameters->time_step);
    if(ts_counter){
      *u_time_average /= ts_counter;
      mpi_driver->Write_Global_Array_To_Disk("u_time_average", 
					*u_time_average, parameters->time_step);
      u_time_average->Set_All_Elements_To(VECTOR_3D<T>(0.,0.,0.));
      ts_counter = 0;
    }
  }

 private:
  CURVILINEAR_GRID<T> *grid;
  PARAMETERS<T> *parameters;
  MPI_DRIVER<T> *mpi_driver;
  ARRAY_3D<VECTOR_3D<T> > *u;
 
  int array_size, host_cpu, ts_counter;
  VECTOR_3D<int> center;
  T *u_center_aggregate;
  ARRAY_3D<VECTOR_3D<T> > *u_time_average;
};
#endif
