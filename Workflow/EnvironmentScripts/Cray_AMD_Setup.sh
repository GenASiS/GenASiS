
module load cpe/25.09
module load PrgEnv-amd

export GENASIS_MACHINE=Cray_AMD

#export OLCF_AFAR_ROOT=/autofs/nccs-svm1_sw/crusher/ums/compilers/afar/rocm-afar-4106
#export OLCF_AFAR_ROOT=/autofs/nccs-svm1_sw/crusher/ums/compilers/afar/rocm-afar-6356-drop-4.1.0
#export OLCF_AFAR_ROOT=/autofs/nccs-svm1_sw/crusher/ums/compilers/afar/rocm-afar-7450-drop-6.0.0
#export OLCF_AFAR_ROOT=/autofs/nccs-svm1_sw/crusher/ums/compilers/afar/rocm-afar-7702-drop-6.1.0
#export OLCF_AFAR_ROOT=/lustre/orion/world-shared/stf006/reubendb/sw/frontier/afar/rocm-afar-8473-drop-7.1.0

AFAR_UMS_BASEDIR=/sw/crusher/ums/compilers/afar
AFAR_UMS_LATEST=$(ls -d --color=never ${AFAR_UMS_BASEDIR}/*/ | tail -n1)
export OLCF_AFAR_ROOT=${AFAR_UMS_LATEST}

export PATH=${OLCF_AFAR_ROOT}/lib/llvm/bin:${PATH}
export LD_LIBRARY_PATH=${OLCF_AFAR_ROOT}/lib:${OLCF_AFAR_ROOT}/lib/llvm/lib:${LD_LIBRARY_PATH}

#module use /sw/crusher/ums/compilers/modulefiles
#module load afar/19.0.0-4106.lua


#export CRAY_MPICH_ROOTDIR=/opt/cray/pe/mpich/8.1.28
#export CRAY_MPICH_PREFIX=${CRAY_MPICH_ROOTDIR}/ofi/amd/5.0
#export CRAY_MPICH_GTL_LIB=${CRAY_MPICH_ROOTDIR}/gtl/lib
#export CRAY_MPICH_INC="-I/autofs/nccs-svm1_sw/crusher/ums/compilers/afar/rocm-afar-3804/include/mpich3.4a2"
export CRAY_MPICH_INC="-I${OLCF_AFAR_ROOT}/include/mpich3.4a2"
#export CRAY_MPICH_LIB="-L${CRAY_MPICH_PREFIX}/lib \
#                        ${PE_MPICH_GTL_DIR_amd_gfx908} \
#                        ${CRAY_PMI_POST_LINK_OPTS} \
#                        -lmpifort_amd -lmpi_amd -lmpi -lpmi -lpmi2"

export CRAY_MPICH_LIB="-L${CRAY_MPICH_PREFIX}/lib \
                        ${CRAY_PMI_POST_LINK_OPTS} \
                        -lmpifort_amd -lmpi_amd -lmpi -lpmi -lpmi2"

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"

unset MPICH_GPU_SUPPORT_ENABLED
