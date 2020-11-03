module purge
module load craype-x86-naples
module load craype-network-infiniband
module load slurm shared

module restore PrgEnv-cray
module load craype-accel-amd-gfx908
module unload cray-mvapich2
module load cray-mvapich2_nogpu
module load rocm/3.8.0
module unload cuda10.2/toolkit

export GENASIS_MACHINE=Cray_CCE

HDF5_DIR=/home/users/coe0021/localsw/poplar/hdf5/1.10.12_cce10
SILO_DIR=/home/users/coe0021/localsw/poplar/silo/4.10.2_cce9.1.3
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${HDF5_DIR}/lib:${SILO_DIR}/lib:${LD_LIBRARY_PATH}
