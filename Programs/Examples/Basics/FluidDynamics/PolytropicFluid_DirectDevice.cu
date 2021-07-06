#include "DirectDevice.h"

#ifdef __cplusplus
extern "C"
  { 
#endif
void ComputeConservedPolytropic_C
       ( double *G, double *E, double *N, double *V_1, 
         double *V_2, double *V_3, int nValues );
void ComputePrimitivePolytropic_C
	     ( double *E, double *G, double *N, 
         double *V_1, double *V_2, double *V_3, int nValues );
void ComputeAuxiliaryPolytropic_C
       ( double *P, double *K, double *N, double *E, 
         double *Gamma, int nValues );
void ComputeEigenspeedsPolytropic_C
       ( double *FEP_1, double *FEP_2, double *FEP_3,
         double *FEM_1, double *FEM_2, double *FEM_3,
         double *CS, double *N, double *V_1, double *V_2, 
         double *V_3, double *P, double *Gamma, int nValues );
void ApplyBoundaryConditionsReflectingPolytropic_C
       ( double *E_E, double *Gamma_E, double *E_I, double *Gamma_I,
         int *nB, int *oBE, int *oBI, int *nSizes );
void ComputeRawFluxesPolytropic_C
       ( double *F_D, double *F_S_1, double *F_S_2, double *F_S_3,
         double *F_S_Dim, double *F_G, double *D, double *S_1, double *S_2, 
         double *S_3, double *G, double *P, double *V_Dim, int nValues );
#ifdef __cplusplus
  }
#endif


__global__ void ComputeConservedPolytropicDeviceKernel 
                  ( double *G, double *E, double *N, double *V_1, 
                    double *V_2, double *V_3, int nValues )
  {
  if ( tiD < nValues )
    {
    G [ tiD ] = E [ tiD ] + 0.5 * N [ tiD ] 
                  * ( V_1 [ tiD ] * V_1 [ tiD ]
                        + V_2 [ tiD ] * V_2 [ tiD ] 
                        + V_3 [ tiD ] * V_3 [ tiD ] );
    }
  }


void ComputeConservedPolytropic_C
       ( double *G, double *E, double *N, double *V_1, 
         double *V_2, double *V_3, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );

  hipLaunchKernelGGL
    ( ( ComputeConservedPolytropicDeviceKernel ), 
      grid_Dim, block_Dim, 0, 0, 
      G, E, N, V_1, V_2, V_3, nValues );
  DeviceSynchronize ( );
  }
  

__global__ void ComputePrimitivePolytropicDeviceKernel
                  ( double *E, double *G, double *N, double *V_1, 
                    double *V_2, double *V_3, int nValues )
  {
  double KE;  	
	
  if ( tiD < nValues )
    {
    KE = 0.5 * N [ tiD ] 
           * ( V_1 [ tiD ] * V_1 [ tiD ]
                 + V_2 [ tiD ] * V_2 [ tiD ]
                 + V_3 [ tiD ] * V_3 [ tiD ] );

    E [ tiD ] = G [ tiD ] - KE;
    if ( E [ tiD ] < 0.0 )
      {
      E [ tiD ] = 0.0;
      G [ tiD ] = KE; 
      }
    }
  }


void ComputePrimitivePolytropic_C
	     ( double *E, double *G, double *N, 
         double *V_1, double *V_2, double *V_3, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );

  hipLaunchKernelGGL
    ( ( ComputePrimitivePolytropicDeviceKernel ), grid_Dim, block_Dim, 0, 0,
      E, G, N, V_1, V_2, V_3, nValues );
  DeviceSynchronize ( );
  }


__global__ void ComputeAuxiliaryPolytropicDeviceKernel 
                  ( double *P, double *K, double *N, double *E, 
                    double *Gamma, int nValues )
  {
  if ( tiD < nValues )
    {
    P [ tiD ] = E [ tiD ] * ( Gamma [ tiD ] - 1.0 );
    if ( N [ tiD ] > 0.0 )
      K [ tiD ] = P [ tiD ] / pow ( N [ tiD ], Gamma [ tiD ] );
    else
      K [ tiD ] = 0.0;
    }
  }


void ComputeAuxiliaryPolytropic_C
       ( double *P, double *K, double *N, double *E, 
         double *Gamma, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ));

  hipLaunchKernelGGL 
    ( ( ComputeAuxiliaryPolytropicDeviceKernel ), 
      grid_Dim, block_Dim, 0, 0,
      P, K, N, E, Gamma, nValues ); 
  DeviceSynchronize ( );
  }
  

__global__ void ComputeEigenspeedsPolytropicDeviceKernel
                  ( double *FEP_1, double *FEP_2, double *FEP_3,
                    double *FEM_1, double *FEM_2, double *FEM_3,
                    double *CS, double *N, double *V_1, double *V_2, 
                    double *V_3, double *P,  double *Gamma, int nValues )
  {
  if ( tiD < nValues )
    {
    if ( N [ tiD ] > 0.0 && P [ tiD ] > 0.0 )
      CS [ tiD ] = sqrt ( Gamma [ tiD ] * P [ tiD ] / N [ tiD ] );
    else
      CS [ tiD ] = 0.0;

    FEP_1 [ tiD ] = V_1 [ tiD ] + CS [ tiD ];
    FEP_2 [ tiD ] = V_2 [ tiD ] + CS [ tiD ];
    FEP_3 [ tiD ] = V_3 [ tiD ] + CS [ tiD ];
    FEM_1 [ tiD ] = V_1 [ tiD ] - CS [ tiD ];
    FEM_2 [ tiD ] = V_2 [ tiD ] - CS [ tiD ];
    FEM_3 [ tiD ] = V_3 [ tiD ] - CS [ tiD ];
    }
  
  }
      

void ComputeEigenspeedsPolytropic_C
       ( double *FEP_1, double *FEP_2, double *FEP_3,
         double *FEM_1, double *FEM_2, double *FEM_3,
         double *CS, double *N, double *V_1, double *V_2, 
         double *V_3, double *P, double *Gamma, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );

  hipLaunchKernelGGL 
    ( ( ComputeEigenspeedsPolytropicDeviceKernel ), 
      grid_Dim, block_Dim, 0, 0, 
      FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3,
      CS, N, V_1, V_2, V_3, P, Gamma, nValues ); 
  DeviceSynchronize ( );
  }


__global__ void ApplyBoundaryConditionsReflectingPolytropic
                  ( double *E_E, double *Gamma_E,
                    double *E_I, double *Gamma_I,
                    int *nB, int *oBE, int *oBI, 
                    int iSize, int jSize, int kSize )
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
  
  if ( kV <= nB [ 3 - 1 ]  &&  jV <= nB [ 2 - 1 ]  &&  iV <= nB [ 1 - 1 ] )
    {
    E_E     [ cI_E ]  = E_I     [ cI_I ];
    Gamma_E [ cI_E ]  = Gamma_I [ cI_I ];
    }
  }

     
void ApplyBoundaryConditionsReflectingPolytropic_C
       ( double *E_E, double *Gamma_E, double *E_I, double *Gamma_I,
         int *nB, int *oBE, int *oBI, int *nSizes )
  {
  int nValues = nSizes [ 3 - 1 ] * nSizes [ 2 - 1 ] * nSizes [ 1 - 1 ];
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ) );
  
  int *d_nB, *d_oBE, *d_oBI;
  
  /*
  hipHostMalloc ( &d_nB,  3 * sizeof ( int ) );
  hipHostMalloc ( &d_oBE, 3 * sizeof ( int ) );
  hipHostMalloc ( &d_oBI, 3 * sizeof ( int ) );
  */
  
  hipMalloc ( &d_nB,  3 * sizeof ( int ) );
  hipMalloc ( &d_oBE, 3 * sizeof ( int ) );
  hipMalloc ( &d_oBI, 3 * sizeof ( int ) );
  
  hipMemcpy ( d_nB, nB, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_oBE, oBE, 3 * sizeof ( int ), hipMemcpyDefault );
  hipMemcpy ( d_oBI, oBI, 3 * sizeof ( int ), hipMemcpyDefault );
  
  hipLaunchKernelGGL
    ( ( ApplyBoundaryConditionsReflectingPolytropic ), 
      grid_Dim, block_Dim, 0, 0,
      E_E, Gamma_E, E_I, Gamma_I, d_nB, d_oBE, d_oBI, 
      nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] );
  
  DeviceSynchronize (  );
  
  hipFree ( d_oBI );
  hipFree ( d_oBE );
  hipFree ( d_nB );
  }
 

__global__ void ComputeRawFluxesPolytropicDeviceKernel
                  ( double *F_D, double *F_S_1, double *F_S_2, 
                    double *F_S_3, double *F_S_Dim, double *F_G,
                    double *D, double *S_1, double *S_2, double *S_3,
                    double *G, double *P, double *V_Dim, int nValues )
  {
  if ( tiD < nValues )
    {
    F_D [ tiD ]     = D [ tiD ]   * V_Dim [ tiD ];
    F_S_1 [ tiD ]   = S_1 [ tiD ] * V_Dim [ tiD ];
    F_S_2 [ tiD ]   = S_2 [ tiD ] * V_Dim [ tiD ];
    F_S_3 [ tiD ]   = S_3 [ tiD ] * V_Dim [ tiD ];
    F_S_Dim [ tiD ] = F_S_Dim [ tiD ] + P [ tiD ];
    F_G [ tiD ]     = ( G [ tiD ] + P [ tiD ] ) * V_Dim [ tiD ];
    } 
  }


void ComputeRawFluxesPolytropic_C
       ( double *F_D, double *F_S_1, double *F_S_2, double *F_S_3,
         double *F_S_Dim, double *F_G, double *D, double *S_1, double *S_2, 
         double *S_3, double *G, double *P, double *V_Dim, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x -1 ) / block_Dim.x ));

  hipLaunchKernelGGL
    ( ( ComputeRawFluxesPolytropicDeviceKernel ), grid_Dim, block_Dim, 0, 0, 
      F_D, F_S_1, F_S_2, F_S_3, F_S_Dim, F_G, D, S_1, S_2, S_3, G, P, 
      V_Dim, nValues ); 
  DeviceSynchronize ( );
  }
