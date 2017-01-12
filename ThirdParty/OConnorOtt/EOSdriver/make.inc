#include $(PATH_GENASIS)/Build/Machines/Makefile_$(GENASIS_MACHINE)
F90=$(FORTRAN_COMPILE) $(FFLAGS) $(FORTRAN_BOUND_CHECK)
HDF5INCS=$(INCLUDE_HDF5)
HDF5LIBS=$(LIBRARY_HDF5)
#-L/home/cott/opt/hdf-1.8.1-gcc4.3/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz
