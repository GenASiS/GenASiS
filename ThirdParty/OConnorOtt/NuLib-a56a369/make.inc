include ../../../../Build/Machines/Makefile_$(GENASIS_MACHINE)
FFLAGS = $(FORTRAN_$(PURPOSE))
F90 := $(shell echo $(FORTRAN_COMPILE) $(FFLAGS) $(FORTRAN_BOUND_CHECK )| sed 's/-c //g')
F77 := $(shell echo $(FORTRAN_COMPILE) $(FORTRAN_FIXEDFORM) $(FFLAGS) $(FORTRAN_BOUND_CHECK )| sed 's/-c //g')
#F90= mpif90
#F90FLAGS= -O3
#F90FLAGS= -g -Og -ffpe-trap=invalid,zero,overflow -fbounds-check

#if you want to take advantage of openmp, set OPENMP=1
OPENMP=0

#if you want to take advantage of mpi, set MPI=1
#mpif90 must be in your PATH (openmpi)
MPI = 0

MODINC="-I ./"

#You must have a HDF5 version installed with
#the _same_ compiler, this often means compiling
#it by yourself.
#See README file for compiling HDF5

HDF5INCS=$(INCLUDE_HDF5)
HDF5LIBS=$(LIBRARY_HDF5)
#HDF5DIR=/usr/local/Cellar/hdf5/1.8.16_1
#HDF5DIR=/usr/local/hdf5-1.8.17
#HDF5INCS=-I$(HDF5DIR)/include
#HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz

#additional flag needed to include other's work, like Hempel's EOS or Timmes' Helmholtz EOS
NUCLEI_HEMPEL = 0
HELMHOLTZ_EOS = 0
WEAK_RATES = 0
