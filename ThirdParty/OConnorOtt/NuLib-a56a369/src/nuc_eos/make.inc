include ../../../../../Build/Machines/Makefile_$(GENASIS_MACHINE)
F90=$(FORTRAN_COMPILE) $(FFLAGS) $(FORTRAN_BOUND_CHECK)
HDF5INCS=$(INCLUDE_HDF5)
HDF5LIBS=$(LIBRARY_HDF5)