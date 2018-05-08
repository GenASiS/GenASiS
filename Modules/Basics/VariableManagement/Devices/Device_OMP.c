#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void * AllocateTargetDouble_OMP ( int nValues )
  {
  int iDevice;
  void * D_Pointer;
  
  iDevice = omp_get_default_device();
  
  /*
  printf("nValues Alloc: %d\n", nValues );
  printf("iDevice: %d\n", iDevice );
  */ 
  D_Pointer = omp_target_alloc ( sizeof ( double ) * nValues, iDevice );
  
  // printf("D_Pointer : %p\n", D_Pointer);
  
  return D_Pointer;
  }


int AssociateTargetDouble_OMP 
      ( void * Host, void * Device, int nValues, int oValue )
  {    
  int iDevice, retval;
  size_t Size, Offset;
  
  /*
  printf("nValues Assoc: %d\n", nValues );
  printf("Host     : %p\n", Host );
  printf("Device   : %p\n", Device);
  */
  
  Size    = sizeof ( double ) * nValues;
  Offset  = sizeof ( double ) * oValue;
  iDevice = omp_get_default_device();
  
  /*
  printf("Host   : %p\n", Host );
  printf("Device : %p\n", Device);
  printf("Size   : %d\n", Size);
  printf("Offset : %d\n", Offset);
  printf("Device : %d\n", iDevice);
  */
  
  retval = omp_target_associate_ptr ( Host, Device, Size, Offset, iDevice );
  //printf("retval assoc: %d\n", retval);
  return retval;
  }
  

void FreeTarget_OMP ( void * D_Pointer )
  {
  int iDevice;
  
  iDevice = omp_get_default_device();
  omp_target_free ( D_Pointer, iDevice );
  }

  
int DisassociateTarget_OMP ( void * Host )
  {
  int iDevice, retval;
  
  iDevice = omp_get_default_device();
  
  retval = omp_target_disassociate_ptr ( Host, iDevice );
  //printf("retval disassoc : %d\n", retval);
  return retval;
  }
  