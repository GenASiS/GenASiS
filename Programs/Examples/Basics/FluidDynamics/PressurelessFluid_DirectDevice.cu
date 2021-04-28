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
void ApplyBoundaryConditionsReflectingPressureless_C
       ( double *N_E, double *VI_E, double *VJ_E, double *VK_E, 
	       double *N_I, double *VI_I, double *VJ_I, double *VK_I,
         int *nB, int *oBE, int *oBI, int *nSizes );
void ComputeRiemannSolverInputPressureless_C
       ( double *AP_I, double *AP_O, double *AM_I, double *AM_O,
         double *LP_I, double *LP_O, double *LM_I, double *LM_O,
         int *lV, int *uV, int *iaS_M, int *iaS_P, int *nSizes );
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


__global__ void ApplyBoundaryConditionsReflectingPressureless
                  ( double *N_E, double *VI_E, double *VJ_E, double *VK_E, 
                    double *N_I, double *VI_I, double *VJ_I, double *VK_I,
                    int *nB, int *oBE, int *oBI, int iSize, int jSize, 
                    int kSize )
  {
  int iV, jV, kV, iI_I, jI_I, kI_I, iI_E, jI_E, kI_E;  // base-1 indexing
  int cI_I, cI_E; // 1D index for Interior and Exterior, base-0 indexing
 
  iV = ( tiD ) % nB [ 1 - 1 ] + 1;
  jV = ( tiD /  nB [ 1 - 1 ] ) % nB [ 2 - 1 ] + 1; 
  kV = tiD / ( nB [ 1 - 1 ] * nB [ 2 - 1 ] ) + 1;
  
  iI_I = oBI [ 0 ] + iV;
  jI_I = oBI [ 1 ] + jV;
  kI_I = oBI [ 2 ] + kV;
	
  cI_I = ( iI_I  +  iSize * ( jI_I - 1 )  
                 +  iSize * jSize * ( kI_I - 1 ) ) - 1;
	
  iI_E = oBE [ 0 ] + iV;
  jI_E = oBE [ 1 ] + jV;
  kI_E = oBE [ 2 ] + kV;
	
  cI_E = ( iI_E  +  iSize * ( jI_E - 1 )  
	               +  iSize * jSize * ( kI_E - 1 ) ) - 1;

  if ( kV <= nB [ 3 - 1 ] && jV <= nB [ 2 - 1 ] && iV <= nB [ 1 - 1 ] )
    {
    N_E  [ cI_E ] 	= N_I  [ cI_I ];
    VI_E [ cI_E ] 	= - VI_I [ cI_I ];
    VJ_E [ cI_E ] 	= VJ_I [ cI_I ];
    VK_E [ cI_E ] 	= VK_I [ cI_I ];
    }
  }


void ApplyBoundaryConditionsReflectingPressureless_C
       ( double *N_E, double *VI_E, double *VJ_E, double *VK_E, 
	       double *N_I, double *VI_I, double *VJ_I, double *VK_I,
         int *nB, int *oBE, int *oBI, int *nSizes )
  {
  int nValues = nSizes [ 3 - 1 ] * nSizes [ 2 - 1 ] * nSizes [ 1 - 1 ];
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );
  
  int *d_nB, *d_oBE, *d_oBI;
  
  hipHostMalloc ( &d_nB,  3 * sizeof ( int ) );
  hipHostMalloc ( &d_oBE, 3 * sizeof ( int ) );
  hipHostMalloc ( &d_oBI, 3 * sizeof ( int ) );
  
  hipMemcpy ( d_nB, nB, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_oBE, oBE, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_oBI, oBI, 3 * sizeof ( int ), hipMemcpyDefault );
  
  hipLaunchKernelGGL 
    ( ( ApplyBoundaryConditionsReflectingPressureless ),
      dim3 ( grid_Dim ), dim3 ( block_Dim ), 0, 0, 
	    N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, d_nB, d_oBE, d_oBI, 
	    nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] ); 

  DeviceSynchronize (  ); 
  
  hipFree ( d_oBI );
  hipFree ( d_oBE );
  hipFree ( d_nB );
  }
  

__global__ void ComputeRiemannSolverInputPressurelessDeviceKernel
                  ( double *AP_I, double *AP_O, double *AM_I, double *AM_O,
                    double *LP_I, double *LP_O, double *LM_I, double *LM_O, 
                    int *lV, int *uV, int *iaS_M, int *iaS_P,
                    int iSize, int jSize, int kSize )
  {
  int iV, jV, kV, iaVS_M, jaVS_M, kaVS_M, iaVS_P, jaVS_P, kaVS_P;
  int A_iD, VS_MiD, VS_PiD; //1D indices
    
  iV = ( tiD ) % ( uV [ 0 ] ) + 1;
  jV = ( tiD / uV [ 0 ] ) %  uV [ 1 ] + 1;
  kV = tiD / ( uV [ 0 ] *  uV [ 1 ] ) + 1;
   
  A_iD = ( iV + iSize * ( jV - 1 )
              + iSize * jSize * ( kV - 1 ) ) - 1;
  
  if ( kV >= lV [ 2 ] && kV <= uV [ 2 ] )
    {
    if ( jV >= lV [ 1 ] && jV <= uV [ 1 ] )
      {
      if ( iV >= lV [ 0 ] && iV <= uV [ 0 ] )
        { 
        iaVS_M = iV + iaS_M [ 0 ];
        jaVS_M = jV + iaS_M [ 1 ];
        kaVS_M = kV + iaS_M [ 2 ];

        VS_MiD = ( iaVS_M + iSize * ( jaVS_M - 1 )
                          + iSize * jSize * ( kaVS_M - 1 ) ) - 1;
        
        iaVS_P = iV + iaS_P [ 0 ];
        jaVS_P = jV + iaS_P [ 1 ];
        kaVS_P = kV + iaS_P [ 2 ];
  
        VS_PiD = ( iaVS_P + iSize * ( jaVS_P - 1 )
                          + iSize * jSize * ( kaVS_P - 1 ) ) - 1; 
  
        AP_I [ A_iD ] 
          = max ( 0.0, max ( + LP_O [ VS_MiD ], + LP_I [ A_iD ] ) );
        AP_O [ A_iD ] 
          = max ( 0.0, max ( + LP_O [ A_iD ], + LP_I [ VS_PiD ] ) );
        AM_I [ A_iD ] 
          = max ( 0.0, max ( - LM_O [ VS_MiD ], - LM_I [ A_iD ] ) );
        AM_O [ A_iD ] 
          = max ( 0.0, max ( - LM_O [ A_iD ], - LM_I [ VS_PiD ] ) );
        }
      }
    }    
  }

      
void ComputeRiemannSolverInputPressureless_C
       ( double *AP_I, double *AP_O, double *AM_I, double *AM_O,
         double *LP_I, double *LP_O, double *LM_I, double *LM_O,
         int *lV, int *uV, int *iaS_M, int *iaS_P, int *nSizes )
  {
  int nValues = nSizes [ 3 - 1 ] * nSizes [ 2 - 1 ] * nSizes [ 1 - 1 ];
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ));
  
  int *d_lV, *d_uV, *d_iaS_M, *d_iaS_P; 
  
  hipHostMalloc ( &d_lV, 3 * sizeof ( int ) );
  hipHostMalloc ( &d_uV, 3 * sizeof ( int ) );
  hipHostMalloc ( &d_iaS_M, 3 * sizeof ( int ) );
  hipHostMalloc ( &d_iaS_P, 3 * sizeof ( int ) );
  
  hipMemcpy ( d_lV, lV, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_uV, uV, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_iaS_M, iaS_M, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_iaS_P, iaS_P, 3 * sizeof ( int ), hipMemcpyDefault );

  hipLaunchKernelGGL
    ( ( ComputeRiemannSolverInputPressurelessDeviceKernel ), 
      grid_Dim, block_Dim, 0, 0,
      AP_I, AP_O, AM_I, AM_O, LP_I, LP_O, LM_I, LM_O, d_lV, d_uV, 
      d_iaS_M, d_iaS_P, nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] );
  DeviceSynchronize ( );
  
  hipFree ( d_iaS_P );
  hipFree ( d_iaS_M );
  hipFree ( d_uV );
  hipFree ( d_lV );
  }
