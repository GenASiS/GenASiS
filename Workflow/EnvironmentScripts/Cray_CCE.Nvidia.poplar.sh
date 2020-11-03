module purge
module load craype-x86-naples
module load PrgEnv-cray
module load craype-network-infiniband
module load slurm
module load craype-accel-nvidia70
module load gcc

HDF5_DIR=/home/users/coe0021/localsw/poplar/hdf5/1.10.12_cce10
SILO_DIR=/home/users/coe0021/localsw/poplar/silo/4.10.2_cce9.1.3
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${HDF5_DIR}/lib:${SILO_DIR}/lib:${LD_LIBRARY_PATH}
