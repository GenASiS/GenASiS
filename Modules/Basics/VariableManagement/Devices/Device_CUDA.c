#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

void * AllocateHostDouble_CUDA ( int nValues )
  {
  void * Host;
  cudaError_t Status;
  
  Status = cudaHostAlloc ( &Host, sizeof ( double ) * nValues, 
                           cudaHostAllocPortable );
  if ( Status != cudaSuccess )
    return NULL;
  
  return Host;
  }


void FreeHost_CUDA ( void * Host )
  {
  cudaError_t Status;
  
  Status = cudaFreeHost ( Host );
  }
