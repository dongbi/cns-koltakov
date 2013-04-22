OBJS = array_3d.o curvilinear_grid.o navier_stokes_solver.o
CC = mpicc
CFLAGS = -c -O3  

test: $(OBJS) test.cpp
	$(CC) $(OBJS) test.cpp -o test

navier_stokes_solver.o : navier_stokes_solver.h navier_stokes_solver.cpp mpi_driver.h parameters.h curvilinear_grid.h curvilinear_moving_grid.h convection.h pressure.h scalar.h potential_energy.h data_aggregator.h tridiagonal_solver.h metric_quantities.h universal_limiter.h moving_grid_engine.h interpolant.h array_3d.h array_2d.h vector_3d.h parameter_file_parser.h
	$(CC) $(CFLAGS) navier_stokes_solver.cpp

curvilinear_grid.o : curvilinear_grid.h curvilinear_grid.cpp array_3d.h array_1d.h parameters.h mpi_driver.h vector_3d.h
	$(CC) $(CFLAGS) curvilinear_grid.cpp

array_3d.o : array_3d.h array_3d.cpp array_2d.h vector_3d.h
	$(CC) $(CFLAGS) array_3d.cpp

run1:
	mpiexec -n 1 ./test
run4:
	mpiexec -n 4 ./test
run8:
	mpiexec -n 8 ./test
run12:
	mpiexec -n 12 ./test
run16:
	mpiexec -n 16 ./test

clean:
	\rm -f *.o *~ test

clobber:
	\rm -rf *.o *~ test ./output/*
