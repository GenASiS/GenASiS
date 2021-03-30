#!/bin/bash

module load gcc/10.2.0
module load cuda hdf5 netlib-lapack

GCC_UMS_DIR=/sw/summit/ums/stf010/gcc
latest=$(ls --color=never ${GCC_UMS_DIR} | tail -n1)
export GCC_ROOT=$GCC_UMS_DIR/$latest

echo "Using GCC in $GCC_ROOT"

export GENASIS_MACHINE=POWER_GCC

export OMPI_CC=${GCC_ROOT}/bin/gcc
export OMPI_CXX=${GCC_ROOT}/bin/g++
export OMPI_FC=${GCC_ROOT}/bin/gfortran
export LD_LIBRARY_PATH=${GCC_ROOT}/lib64:${LD_LIBRARY_PATH}

echo "mpif90 --version is: "
mpif90 --version
