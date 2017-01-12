#!/bin/bash
set +x
export GENASIS_MACHINE=$1

echo "Unloading modules to start clean ... "

if [[ $GENASIS_MACHINE == Titan* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale
fi

if [[ $GENASIS_MACHINE == Jaguar* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale xt-asyncpe
fi

if [[ $GENASIS_MACHINE == Chester* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale
fi

if [[ $GENASIS_MACHINE == Hopper* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale
fi

if [[ $GENASIS_MACHINE == Interlagos* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale
fi

if [[ $GENASIS_MACHINE == Smoky* ]]; then
  module unload PE-pgi PE-gnu PE-intel PrgEnv-pathscale
  module unload pgi gcc cce pathscale petsc
  module unload fftw hdf5 petsc silo
fi

if [[ $GENASIS_MACHINE == Kraken* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale
fi

if [[ $GENASIS_MACHINE == Mars* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray
fi

if [[ $GENASIS_MACHINE == Metis* ]]; then
  module unload fftw hdf5 petsc silo
  module unload pgi gcc cce pathscale petsc
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray
fi

if [[ $GENASIS_MACHINE == Darter* ]]; then
  module unload fftw cray-hdf5 cray-petsc
  module unload pgi gcc cce intel
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-intel
fi

if [[ $GENASIS_MACHINE == Comet* ]]; then
  module purge
fi

#if [[ $GENASIS_MACHINE == Alpha* ]]; then
#fi


echo "Loading environment for $GENASIS_MACHINE"

#-- titan.ccs.ornl.gov 
if [[ "$GENASIS_MACHINE" == "Titan_Cray" ]]; then
  module load PrgEnv-cray
  module load cray-petsc
  module load fftw
  module load ddt
fi

#-- jaguar.ccs.ornl.gov 
if [[ "$GENASIS_MACHINE" == "Jaguar_Cray" ]]; then
  module load PrgEnv-cray
  module swap cce cce/8.1.0.160
  module load fftw
  module load petsc
  module load ddt/3.1-22442
fi

if [[ "$GENASIS_MACHINE" == "Jaguar_PGI" ]]; then
  module load PrgEnv-pgi
  module swap pgi pgi/12.2.0
  module load fftw
  module load ddt
fi

if [[ "$GENASIS_MACHINE" == "Chester_Cray" ]]; then
  module load PrgEnv-cray
  module swap cce cce/8.1.2
  module load petsc
  module load fftw
  module load ddt
fi

if [[ "$GENASIS_MACHINE" == "Chester_PGI" ]]; then
  module load PrgEnv-pgi
  module load petsc
  module load fftw
  module load ddt
fi

if [[ "$GENASIS_MACHINE" == "Interlagos_Cray" ]]; then
  module load subversion/1.6.17
  module load PrgEnv-cray
  module unload cce
  module load cce/8.0.0.135
  module load petsc
  module load fftw
  module load ddt
fi

if [[ "$GENASIS_MACHINE" == "Interlagos_PGI" ]]; then
  module load subversion/1.6.17
  module load PrgEnv-pgi
  module load petsc
  module load fftw
  module load ddt
fi

if [[ "$GENASIS_MACHINE" == "Mars_Cray" ]]; then
  module load PrgEnv-cray
  module load cray-hdf5
  module load cray-petsc
  module load fftw
  module load silo
fi

if [[ "$GENASIS_MACHINE" == "Mars_GNU" ]]; then
  module load PrgEnv-gnu
  module load cray-hdf5
  module load cray-petsc
  module load fftw
  module load silo
fi

if [[ "$GENASIS_MACHINE" == "Metis_Cray" ]]; then
  module load PrgEnv-cray
  module load cray-hdf5
  module load fftw
  module load silo/4.10.2
fi

if [[ "$GENASIS_MACHINE" == "Metis_PGI" ]]; then
  module load PrgEnv-pgi
  module load cray-hdf5
  module load fftw
  module load silo/4.10.2
fi

if [[ "$GENASIS_MACHINE" == "Metis_GNU" ]]; then
  module load PrgEnv-gnu
  module swap gcc gcc/6.1.0
  module load cray-hdf5
  module load fftw
  module load silo/4.10.2
fi

if [[ "$GENASIS_MACHINE" == "Darter_Cray" ]]; then
  module load PrgEnv-cray
  module load cray-petsc
  module load fftw
  module load silo/4.10.2
  module load python/.2.7.11
  export PYTHONWARNINGS="ignore:Unverified HTTPS request"
fi

if [[ "$GENASIS_MACHINE" == "Darter_GNU" ]]; then
  module load PrgEnv-gnu
  module swap gcc gcc/5.3.0
  module load cray-petsc
  module load fftw
  module load silo/4.10.2
  export PYTHONWARNINGS="ignore:Unverified HTTPS request"
fi

if [[ "$GENASIS_MACHINE" == "Darter_Intel" ]]; then
  module load PrgEnv-intel
  module swap intel /nics/d/home/rbudiard/software/modules/intel/14.0.1.106.fixed
  module load cray-hdf5
  module load cray-petsc
  module load fftw
  module load silo
  module load python/.2.7.11
  export PYTHONWARNINGS="ignore:Unverified HTTPS request"
fi

if [[ "$GENASIS_MACHINE" == "Beacon_Intel" ]]; then
  source ${MODULESHOME}/init/bash
  module load silo/4.10
  module load hdf5
fi

if [[ "$GENASIS_MACHINE" == "Beacon_GNU" ]]; then
  source ${MODULESHOME}/init/bash
  module purge
  module unuse /sw/beacon/modulefiles
  module use   /sw/cs400_centos7.2_pe2016-08/modulefiles
  module load PE-gnu silo hdf5
  module load /global/opt/modulefiles/torque/4.2.6
fi

#-- hopper.nersc.gov
if [[ "$GENASIS_MACHINE" == "Hopper_Cray" ]]; then
  module load PrgEnv-cray
  module swap cce cce/7.2.7
  module load fftw
fi

if [[ "$GENASIS_MACHINE" == "Smoky_PGI" ]]; then
  module load PE-pgi
  module swap pgi pgi/11.10
  module load blas-goto
  module load petsc
  module load fftw
  module load ddt
  module swap subversion subversion/1.6.12
fi

if [[ "$GENASIS_MACHINE" == "Smoky_Intel" ]]; then
  module load PE-intel
  module unload intel
  module load /lustre/widow2/proj/ast005/software/modulefiles/intel/12
  module load petsc
  module load ddt
  module swap subversion subversion/1.6.12
fi

if [[ "$GENASIS_MACHINE" == "Comet_Intel" ]]; then
  export MODULEPATH=/oasis/projects/nsf/nic103/rbudiard/software/swtree/comet/modulefiles:$MODULEPATH
  export PE_ADD_RPATH=1
  module load PE-intel
  module load hdf5/1.8.17
  module load silo
fi

if [[ "$GENASIS_MACHINE" == "CCondo_GNU" ]]; then
  export PE_ADD_RPATH=1
  module load PE-gnu
  module swap gcc gcc/6.2.0
  module swap openmpi
  module load silo/4.10.2
fi

#-- Christian's laptop
if [[ "$GENASIS_MACHINE" == "Alpha_Intel" ]]; then
#  source /Users/cca/opt/intel/composer_xe_2011_sp1.9.289/bin/compilervars.sh intel64
  source /Users/cca/opt/intel/composer_xe_2013.1.119/bin/compilervars.sh intel64
fi
