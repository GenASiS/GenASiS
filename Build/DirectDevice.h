#include <stdio.h>
#include <math.h> 

#ifdef __HIP_PLATFORM_HCC__
#include <hip/hip_runtime.h> 
#define DeviceSynchronize hipDeviceSynchronize
#else
#include <cuda_runtime_api.h>
#include <cuda.h>
#define DeviceSynchronize cudaDeviceSynchronize
#endif

#define BLOCK_DIM 16
#define tiD ( blockIdx.x * blockDim.x + threadIdx.x )

