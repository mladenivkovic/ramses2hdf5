MODOBJ= 

F90=gfortran
FFLAGS=-ffree-line-length-none -Ofast -Wall -x f95-cpp-input -Wall -fbacktrace -g #  -std=f2008
LIBFLAGS= -I$(HDF5_ROOT)/include -L$(HDF5_ROOT)/lib -lhdf5_fortran

all: ramses2hdf5

ramses2hdf5: ramses2hdf5.f90 $(MODOBJ)
	$(F90) $(FFLAGS) $^ -o $@ $(LIBFLAGS)




%.o: %.f90
	$(F90) $(FFLAGS) -c $^ -o $@
clean:
	rm *.o *.mod ramses2hdf5
