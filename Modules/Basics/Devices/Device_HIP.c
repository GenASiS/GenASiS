#include <stdlib.h>
#include <stdio.h>

#ifdef ENABLE_OMP_OFFLOAD
#include <hip/hip_runtime_api.h>
#endif

void * AllocateHostDouble_Device ( int nValues )
  {
  void * Host;

  #ifdef ENABLE_OMP_OFFLOAD
  hipError_t Status;
  int Flag;
  
  Flag = 0;
  Status = hipHostMalloc ( &Host, sizeof ( double ) * nValues, 
                           Flag );
  if ( Status != hipSuccess )
    return NULL;
  #else
  Host = malloc ( sizeof ( double ) * nValues );
  #endif
  return Host;
  }


void FreeHost_Device ( void * Host )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  hipError_t Status;
  
  Status = hipHostFree ( Host );
  #else
  free ( Host );
  #endif
  }
  

int DeviceMemGetInfo_Device ( size_t * Free, size_t * Total )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  
  hipError_t Status;
  Status = hipMemGetInfo ( Free, Total );
  
  if ( Status != hipSuccess )
    return (int) Status;
  else
    return 0;
  #else
  Free  = 0;
  Total = 0;
  return -1;
  #endif
  }
