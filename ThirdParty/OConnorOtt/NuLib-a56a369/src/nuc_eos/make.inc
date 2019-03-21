include ../../../../../Build/Machines/Makefile_$(GENASIS_MACHINE)
F90 := $(shell echo $(FORTRAN_COMPILE) $(FFLAGS) $(FORTRAN_BOUND_CHECK )| sed 's/-c //g')
HDF5INCS=$(INCLUDE_HDF5)
HDF5LIBS=$(LIBRARY_HDF5)
