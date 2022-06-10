#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <omp.h>

int OnTarget_OMP ( void * Host )
  {
  int iDevice;
  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device();
  return omp_target_is_present ( Host, iDevice );
  #else
  return 0;
  #endif
  }


void * AllocateTargetInteger_OMP ( int nValues )
  {
  int iDevice;
  void * D_Pointer;
  
  D_Pointer = NULL; 
  
  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device();
  
  /*
  printf("nValues Alloc: %d\n", nValues );
  printf("iDevice: %d\n", iDevice );
  */ 
  D_Pointer = omp_target_alloc ( sizeof ( int ) * nValues, iDevice );
  
  // printf("D_Pointer : %p\n", D_Pointer);
  #endif 
  
  return D_Pointer;
  }


void * AllocateTargetDouble_OMP ( int nValues )
  {
  int iDevice;
  void * D_Pointer;
  
  D_Pointer = NULL; 
  //printf("allocate target double\n");
  
  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device();
  
  /*
  printf("nValues Alloc: %d\n", nValues );
  printf("pre iDevice: %d\n", iDevice );
  */ 
  D_Pointer = omp_target_alloc ( sizeof ( double ) * nValues, iDevice );
  
  //printf("D_Pointer : %p\n", D_Pointer);
  
  //omp_set_default_device(iDevice);
  /*
  iDevice = omp_get_default_device();
  printf("post iDevice: %d\n", iDevice );
  */
  // printf("D_Pointer : %p\n", D_Pointer);
  #endif 
  
  return D_Pointer;
  }


void * AllocateTargetLogical_OMP ( int nValues )
  {
  int iDevice;
  void * D_Pointer;
  
  D_Pointer = NULL; 
  
  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device();
  
  /*
  printf("nValues Alloc: %d\n", nValues );
  printf("iDevice: %d\n", iDevice );
  */ 
  D_Pointer = omp_target_alloc ( sizeof ( _Bool ) * nValues, iDevice );
  
  // printf("D_Pointer : %p\n", D_Pointer);
  #endif 
  
  return D_Pointer;
  }


int AssociateTargetInteger_OMP 
      ( void * Host, void * Device, int nValues, int oValue )
  {    
  int iDevice, retval;
  size_t Size, Offset;
  
  retval = -1;
  
  //printf("nValues Assoc: %d\n", nValues );

  #ifdef ENABLE_OMP_OFFLOAD
  Size    = sizeof ( int ) * nValues;
  Offset  = sizeof ( int ) * oValue;
  iDevice = omp_get_default_device();
  
  /*
  printf("Host   : %p\n", Host );
  printf("Device : %p\n", Device);
  printf("Size   : %d\n", Size);
  printf("Offset : %d\n", Offset);
  printf("Device pre : %d\n", iDevice);
  */
  retval = omp_target_associate_ptr ( Host, Device, Size, Offset, iDevice );
  
  //iDevice = omp_get_default_device();
  //printf("Device post : %d\n", iDevice);
  
  #endif 
  //printf("retval assoc: %d\n", retval);
  return retval;
  }
  

int AssociateTargetDouble_OMP 
      ( void * Host, void * Device, int nValues, int oValue )
  {    
  int iDevice, retval;
  size_t Size, Offset;
  
  retval = -1;
  
  //printf("nValues Assoc: %d\n", nValues );

  #ifdef ENABLE_OMP_OFFLOAD
  Size    = sizeof ( double ) * nValues;
  Offset  = sizeof ( double ) * oValue;
  iDevice = omp_get_default_device();
  
  /*
  printf("Host   : %p\n", Host );
  printf("Device : %p\n", Device);
  printf("Size   : %d\n", Size);
  printf("Offset : %d\n", Offset);
  printf("Device pre : %d\n", iDevice);
  */
  retval = omp_target_associate_ptr ( Host, Device, Size, Offset, iDevice );
  
  //iDevice = omp_get_default_device();
  //printf("Device post : %d\n", iDevice);
  
  #endif 
  //printf("retval assoc: %d\n", retval);
  return retval;
  }
  

int AssociateTargetLogical_OMP 
      ( void * Host, void * Device, int nValues, int oValue )
  {    
  int iDevice, retval;
  size_t Size, Offset;
  
  retval = -1;
  
  /*
  printf("nValues Assoc: %d\n", nValues );
  printf("Host     : %p\n", Host );
  printf("Device   : %p\n", Device);
  */

  #ifdef ENABLE_OMP_OFFLOAD
  Size    = sizeof ( _Bool ) * nValues;
  Offset  = sizeof ( _Bool ) * oValue;
  iDevice = omp_get_default_device();
  
  
  printf("Host   : %p\n", Host );
  printf("Device : %p\n", Device);
  printf("Size   : %d\n", Size);
  printf("Offset : %d\n", Offset);
  printf("Device : %d\n", iDevice);
  
  
  retval = omp_target_associate_ptr ( Host, Device, Size, Offset, iDevice );
  #endif 
  //printf("retval assoc: %d\n", retval);
  return retval;
  }
  

void FreeTarget_OMP ( void * D_Pointer )
  {
  int iDevice;
  
  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device();
  omp_target_free ( D_Pointer, iDevice );
  #endif 
  
  }

  
int DisassociateTarget_OMP ( void * Host )
  {
  int iDevice, retval;
  
  retval = -1;

  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device();
  
  retval = omp_target_disassociate_ptr ( Host, iDevice );
  //printf("retval disassoc : %d\n", retval);
  #endif
  
  return retval;
  }


int HostToDeviceCopyDouble_OMP ( void * Host, void * Device, int nValues, 
                             int oHostValue, int oDeviceValue )
  {
  int iHost, iDevice, retval;
  size_t Length, oHV, oDV;
  
  retval = -1;

  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device ( );
  iHost   = omp_get_initial_device ( );
  Length  = sizeof ( double ) * nValues;
  oHV     = sizeof ( double ) * oHostValue;
  oDV     = sizeof ( double ) * oDeviceValue;
  
  //printf("HostToDevice iDevice: %d\n", iDevice);
  
  retval = omp_target_memcpy 
             ( Device, Host, Length, oDV, oHV, iDevice, iHost ); 
  #endif
  
  return retval;
  }

int DeviceToHostCopyDouble_OMP ( void * Device, void * Host, int nValues, 
                             int oDeviceValue, int oHostValue )
  {
  int iHost, iDevice, retval;
  size_t Length, oHV, oDV;
  
  retval = -1;
  
  #ifdef ENABLE_OMP_OFFLOAD
  iDevice = omp_get_default_device ( );
  iHost   = omp_get_initial_device ( );
  Length  = sizeof ( double ) * nValues;
  oHV     = sizeof ( double ) * oHostValue;
  oDV     = sizeof ( double ) * oDeviceValue;
  
  //printf("DeviceToHost iDevice: %d\n", iDevice);
  
  retval = omp_target_memcpy 
             ( Host, Device, Length, oHV, oDV, iHost, iDevice ); 
  #endif
  
  return retval;
  }
  
bool OffloadEnabled ( )
  {
  #ifdef ENABLE_OMP_OFFLOAD
  return true;
  #else
  return false;
  #endif
  }
