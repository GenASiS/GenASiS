#include <stdlib.h>
#include <stdio.h>

#ifdef ENABLE_OMP_OFFLOAD
#include <cuda.h>
#include <cuda_runtime_api.h>
#endif

int SetDevice ( int iDevice )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  cudaError_t Status;
  
  Status = cudaSetDevice ( iDevice );
  if ( Status == cudaSuccess )
    return 0;
  else
    return -1;
  #else
  return -1;
  #endif
  }

int GetDevice ( int * iDevice )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  cudaError_t Status;
  
  Status = cudaGetDevice ( iDevice );
  if ( Status == cudaSuccess )
    return 0;
  else
    return -1;
  #else
  return -1;
  #endif
  }


void * AllocateHostDouble_Device ( int nValues )
  {
  void * Host;

  #ifdef ENABLE_OMP_OFFLOAD
  cudaError_t Status;
  
  Status = cudaHostAlloc ( &Host, sizeof ( double ) * nValues, 
                           cudaHostAllocPortable );
  if ( Status != cudaSuccess )
    return NULL;
  #else
  Host = malloc ( sizeof ( double ) * nValues );
  #endif
  return Host;
  }


void FreeHost_Device ( void * Host )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  cudaError_t Status;
  
  Status = cudaFreeHost ( Host );
  #else
  free ( Host );
  #endif
  }
  

int DeviceMemGetInfo_Device ( size_t * Free, size_t * Total )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  
  cudaError_t Status;
  Status = cudaMemGetInfo ( Free, Total );
  
  if ( Status != cudaSuccess )
    return (int) Status;
  else
    return 0;
  #else
  Free  = 0;
  Total = 0;
  return -1;
  #endif
  }
