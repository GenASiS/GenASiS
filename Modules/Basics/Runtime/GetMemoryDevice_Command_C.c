#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void * Get ( int nValues )
  {
  int iDevice;
  double *D_Pointer;
  
  iDevice = omp_get_default_device();
  D_Pointer = omp_target_alloc ( sizeof ( double ) * nValues, iDevice );
  
  return D_Pointer;
  }
  