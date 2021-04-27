#include "DirectDevice.h"

#ifdef __cplusplus
extern "C"
  { 
#endif
void ComputeConservedPressureless_C
       ( double *D, double *S_1, double *S_2, double *S_3, 
         double *N, double *V_1, double *V_2, double *V_3, int nValues );
void ComputePrimitivePressureless_C
       ( double *N, double *V_1, double *V_2, double *V_3, 
         double *D, double *S_1, double *S_2, double *S_3, int nValues );
void ComputeEigenspeedsPressureless_C
       ( double *FEP_1, double *FEP_2, double *FEP_3,
         double *FEM_1, double *FEM_2, double *FEM_3,
         double *V_1, double *V_2, double *V_3, int nValues );
#ifdef __cplusplus
  }
#endif


__global__ void ComputeConservedPressurelessDeviceKernel 
                  ( double *D, double *S_1, double *S_2, double *S_3, 
                    double *N, double *V_1, double *V_2, double *V_3, int nValues )
  {
  if ( tiD < nValues )
    {	
    D   [ tiD ] = N [ tiD ];
    S_1 [ tiD ] = N [ tiD ] * V_1 [ tiD ];
    S_2 [ tiD ] = N [ tiD ] * V_2 [ tiD ];
    S_3 [ tiD ] = N [ tiD ] * V_3 [ tiD ]; 
    }
  }


void ComputeConservedPressureless_C
       ( double *D, double *S_1, double *S_2, double *S_3, 
         double *N, double *V_1, double *V_2, double *V_3, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );

  hipLaunchKernelGGL 
    ( ( ComputeConservedPressurelessDeviceKernel ), grid_Dim, block_Dim, 0, 0,
      D, S_1, S_2, S_3, N, V_1, V_2, V_3, nValues );
  DeviceSynchronize ( );
  }


__global__ void ComputePrimitivePressurelessDeviceKernel
	                ( double *N, double *V_1, double *V_2, double *V_3, 
	                  double *D, double *S_1, double *S_2, double *S_3, int nValues )
  { 
  if ( tiD < nValues )
    {
    N [ tiD ]  = D [ tiD ];
    if ( N [ tiD ] > 0.0 )
      {
      V_1 [ tiD ] = S_1 [ tiD ] / N [ tiD ];
      V_2 [ tiD ] = S_2 [ tiD ] / N [ tiD ];
      V_3 [ tiD ] = S_3 [ tiD ] / N [ tiD ];
      }
     else
      {
      N   [ tiD ] = 0.0;
      V_1 [ tiD ] = 0.0;
      V_2 [ tiD ] = 0.0;
      V_3 [ tiD ] = 0.0;
      D   [ tiD ] = 0.0;
      S_1 [ tiD ] = 0.0;
      S_2 [ tiD ] = 0.0;
      S_3 [ tiD ] = 0.0;
      }
    }
  }


void ComputePrimitivePressureless_C
       ( double *N, double *V_1, double *V_2, double *V_3, 
         double *D, double *S_1, double *S_2, double *S_3, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );

  hipLaunchKernelGGL
    ( ( ComputePrimitivePressurelessDeviceKernel ), 
      grid_Dim, block_Dim, 0, 0, 
      N, V_1, V_2, V_3, D, S_1, S_2, S_3, nValues );
  DeviceSynchronize ( );
  }


__global__ void ComputeEigenspeedsPressurelessDeviceKernel
	                ( double *FEP_1, double *FEP_2, double *FEP_3,
	                  double *FEM_1, double *FEM_2, double *FEM_3,
	                  double *V_1, double *V_2, double *V_3, int nValues )
  {
  if ( tiD < nValues )
    {
    FEP_1 [ tiD ] = V_1 [ tiD ];
    FEP_2 [ tiD ] = V_2 [ tiD ];
    FEP_3 [ tiD ] = V_3 [ tiD ];
    FEM_1 [ tiD ] = V_1 [ tiD ];
    FEM_2 [ tiD ] = V_2 [ tiD ];
    FEM_3 [ tiD ] = V_3 [ tiD ];
    }
  }


void ComputeEigenspeedsPressureless_C
       ( double *FEP_1, double *FEP_2, double *FEP_3,
         double *FEM_1, double *FEM_2, double *FEM_3,
         double *V_1, double *V_2, double *V_3, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );

  hipLaunchKernelGGL 
    ( ( ComputeEigenspeedsPressurelessDeviceKernel ), 
      grid_Dim, block_Dim, 0, 0, 
      FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, 
      V_1, V_2, V_3, nValues ); 
  DeviceSynchronize ( );
  }

