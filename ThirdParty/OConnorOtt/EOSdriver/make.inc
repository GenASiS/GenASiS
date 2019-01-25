#include $(PATH_GENASIS)/Build/Machines/Makefile_$(GENASIS_MACHINE)
#F90=$(FORTRAN_COMPILE) $(FFLAGS) $(FORTRAN_BOUND_CHECK)
#HDF5INCS=$(INCLUDE_HDF5)
#HDF5LIBS=$(LIBRARY_HDF5)
F90 = mpif90 -g -framework Accelerate -fbacktrace -ffpe-trap=invalid,zero -ffpe-summary=invalid,zero -fbounds-check
HDF5INCS = -I/usr/local/Cellar/hdf5/1.8.16_1/include
HDF5LIBS = -L/usr/local/Cellar/hdf5/1.8.16_1/lib -lhdf5_fortran -lhdf5
#-L/home/cott/opt/hdf-1.8.1-gcc4.3/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz
