#include "DirectDevice.h"

#ifdef __cplusplus
extern "C"
  { 
#endif
  void ComputeDifferences_C
         ( double *V, int *lV, int *uV, int *iaS, int iD, 
           double *dV_Left, double *dV_Right, int *nSizes );
  void ComputeReconstruction_C
       ( double *V, double *dV_Left, double *dV_Right, double Theta,
         double *V_Inner, double *V_Outer, int nValues );
  void ComputeFluxes_C
       ( double *AP_I, double *AP_O, double *AM_I, double *AM_O,
         double *RF_I, double *RF_O, double *U_I, double *U_O,  
         int *lV, int *uV, int *iaS, int iD, double *F_I, double *F_O,
         int *nSizes );
  void ComputeUpdate_C
         ( double *dU, double *F_I, double *F_O,
           double V, double A, double dT, int nValues );
  void AddUpdate_C ( double *O, double *U, double *C, int nValues );
  void CombineUpdates_C ( double *C, double *O, double *U, int nValues );
#ifdef __cplusplus
  }
#endif


__global__ void ComputeDifferencesLeftDeviceKernel
                  ( double *V, int *lV, int *uV, int *iaS, double *dV,
                    int iSize, int jSize, int kSize )
  { 
  int iV, jV, kV, iaVS, jaVS, kaVS;
  int dV_iD, VS_iD; //1D index
     
  iV = ( tiD ) % ( uV [ 0 ] ) + 1;
  jV = ( tiD / uV [ 0 ] ) %  uV [ 1 ] + 1 ;
  kV = tiD / ( uV [ 0 ] *  uV [ 1 ] ) + 1;
     
  dV_iD = ( iV + iSize * ( jV - 1 ) 
               + iSize * jSize * ( kV - 1 ) ) - 1;
  
  if ( kV >= lV [ 2 ] && kV <= uV [ 2 ] )
    {
    if ( jV >= lV [ 1 ] && jV <= uV [ 1 ] )
      {
      if ( iV >= lV [ 0 ] && iV <= uV [ 0 ] )
        { 
        iaVS = iV + iaS [ 0 ];
        jaVS = jV + iaS [ 1 ];
        kaVS = kV + iaS [ 2 ];
  
        VS_iD = ( iaVS + iSize * ( jaVS - 1 )
                       + iSize * jSize * ( kaVS - 1 ) ) - 1;
   
        dV [ dV_iD ] = V [ dV_iD ] - V [ VS_iD ];
        }
      }
    }
  
  }

   
__global__ void ComputeDifferencesRightDeviceKernel
                  ( double *V, int *lV, int *uV, int *iaS, double *dV,
                    int iSize, int jSize, int kSize )
  {
  int iV, jV, kV, iaVS, jaVS, kaVS;
  int dV_iD, VS_iD; //1D index

  
  iV = ( tiD ) % ( uV [ 0 ] ) + 1;
  jV = ( tiD / uV [ 0 ] ) %  uV [ 1 ] + 1;
  kV = tiD / ( uV [ 0 ] *  uV [ 1 ] ) + 1;
   
  dV_iD = ( iV + iSize * ( jV - 1 ) 
               + iSize * jSize * ( kV - 1 ) ) - 1;
  

  if ( kV >= lV [ 2 ] && kV <= uV [ 2 ] )
    {
    if ( jV >= lV [ 1 ] && jV <= uV [ 1 ] )
      {
      if ( iV >= lV [ 0 ] && iV <= uV [ 0 ] )
        { 
        iaVS = iV + iaS [ 0 ];
        jaVS = jV + iaS [ 1 ];
        kaVS = kV + iaS [ 2 ];
    
        VS_iD = ( iaVS + iSize * ( jaVS - 1 )
                       + iSize * jSize * ( kaVS - 1 ) ) - 1;  
  
         dV [ dV_iD ] = V [ VS_iD ] - V [ dV_iD ];
        }
      }
    }
 
  }
 

void ComputeDifferences_C
       ( double *V, int *lV, int *uV, int *iaS, int iD, 
         double *dV_Left, double *dV_Right, int *nSizes )
  {
  int nValues = nSizes [ 3 - 1 ] * nSizes [ 2 - 1 ] * nSizes [ 1 - 1 ];
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x - 1 ) / block_Dim.x ) );
   
  iaS [ iD - 1 ] = -1;
  hipLaunchKernelGGL 
    ( ( ComputeDifferencesLeftDeviceKernel ), 
      dim3 ( grid_Dim ), dim3 ( block_Dim ), 0, 0, 
      V, lV, uV, iaS, dV_Left, nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] );
  
  DeviceSynchronize (  );
  
  iaS [ iD - 1 ] = +1;
  hipLaunchKernelGGL 
    ( ( ComputeDifferencesRightDeviceKernel ), 
        dim3 ( grid_Dim ), dim3 ( block_Dim ), 0, 0,
      V, lV, uV, iaS, dV_Right, nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] );

  DeviceSynchronize (  );
  }
  

__global__ void ComputeReconstructionDeviceKernel
                  ( double *V, double *dV_Left, double *dV_Right, 
                    double Theta, double *V_Inner, double *V_Outer, 
                    int nValues )
  {
  double dV;
  
  if ( tiD < nValues )
    {
    dV = ( copysign ( 0.5, dV_Left [ tiD ] )
             + copysign ( 0.5, dV_Right [ tiD ] ) )
         * min ( min ( abs ( Theta * dV_Left [ tiD ] ),
                       abs ( Theta * dV_Right [ tiD ] ) ),
                 abs ( 0.5 * ( dV_Left [ tiD ] + dV_Right [ tiD ] ) ) );
    V_Inner [ tiD ] = V [ tiD ] - 0.5 * dV;  
    V_Outer [ tiD ] = V [ tiD ] + 0.5 * dV;
    }
  }


void ComputeReconstruction_C
       ( double *V, double *dV_Left, double *dV_Right, double Theta,
         double *V_Inner, double *V_Outer, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x - 1 ) / block_Dim.x ) );
  
  hipLaunchKernelGGL
    ( ( ComputeReconstructionDeviceKernel ), grid_Dim, block_Dim, 0, 0,
        V, dV_Left, dV_Right, Theta, V_Inner, V_Outer, nValues );

  DeviceSynchronize ( );
  }
  
  
__global__ void ComputeFluxesInnerDeviceKernel
                  ( double *AP, double *AM, double *RF_I, double *RF_O, 
                    double *U_I, double *U_O, int *lV, int *uV, int *iaS,
                    double *F, int iSize, int jSize, int kSize )
  {
  int iV, jV, kV, iaVS, jaVS, kaVS;
  int F_iD, VS_iD; //1D index

  iV = ( tiD ) % ( uV [ 0 ] ) + 1;
  jV = ( tiD / uV [ 0 ] ) %  uV [ 1 ] + 1 ;
  kV = tiD / ( uV [ 0 ] *  uV [ 1 ] ) + 1;
   
  F_iD = ( iV + iSize * ( jV - 1 )
              + iSize * jSize * ( kV - 1 ) ) - 1;

  if ( kV >= lV [ 2 ] && kV <= uV [ 2 ] )
    {
    if ( jV >= lV [ 1 ] && jV <= uV [ 1 ] )
      {
      if ( iV >= lV [ 0 ] && iV <= uV [ 0 ] )
        { 
        iaVS = iV + iaS [ 0 ];
        jaVS = jV + iaS [ 1 ];
        kaVS = kV + iaS [ 2 ];

        VS_iD = ( iaVS + iSize * ( jaVS - 1 )
                       + iSize * jSize * ( kaVS - 1 ) ) - 1;
  
        if ( AP [ F_iD ] + AM [ F_iD ] > 0.0 )
          {
          F [ F_iD ] = ( AP [ F_iD ] * RF_O [ VS_iD ]
                            + AM [ F_iD ] * RF_I [ F_iD ]
                            - AP [ F_iD ] * AM [ F_iD ]
                         * ( U_I [ F_iD ]  -  U_O [ VS_iD ] ) )
                       / ( AP [ F_iD ] + AM [ F_iD ] );
          }
        else
          {
          F [ F_iD ] = 0.0;
          }
        }
      }
    } 
  }


__global__ void ComputeFluxesOuterDeviceKernel
                  ( double *AP, double *AM, double *RF_I, double *RF_O,
                    double *U_I, double *U_O, int *lV, int *uV, int *iaS,
                    double *F, int iSize, int jSize, int kSize )
  {
  int iV, jV, kV, iaVS, jaVS, kaVS;
  int F_iD, VS_iD; //1D index
    
    
  iV = ( tiD ) % ( uV [ 0 ] ) + 1;
  jV = ( tiD / uV [ 0 ] ) %  uV [ 1 ] + 1;
  kV = tiD / ( uV [ 0 ] *  uV [ 1 ] ) + 1;
   
  F_iD = ( iV + iSize * ( jV - 1 )
              + iSize * jSize * ( kV - 1 ) ) - 1;
  
  if ( kV >= lV [ 2 ] && kV <= uV [ 2 ] )
    {
    if ( jV >= lV [ 1 ] && jV <= uV [ 1 ] )
      {
      if ( iV >= lV [ 0 ] && iV <= uV [ 0 ] )
        { 
        iaVS = iV + iaS [ 0 ];
        jaVS = jV + iaS [ 1 ];
        kaVS = kV + iaS [ 2 ];
  
        VS_iD = ( iaVS + iSize * ( jaVS - 1 )
                       + iSize * jSize * ( kaVS - 1 ) ) - 1;
  
        if ( AP [ F_iD ] + AM [ F_iD ] > 0.0 )
          {
          F  [ F_iD ] = ( AP [ F_iD ] * RF_O [ F_iD ]
                             + AM [ F_iD ] * RF_I [ VS_iD ]
                             - AP [ F_iD ] * AM [ F_iD ]   
                          * ( U_I [ VS_iD ]  -  U_O [ F_iD ] ) )
                        / ( AP [ F_iD ] + AM [ F_iD ] );
          }
        else
          {
          F [ F_iD ] = 0.0;
          }
        }
      }  
    }    
  
  }      


void ComputeFluxes_C
       ( double *AP_I, double *AP_O, double *AM_I, double *AM_O,
         double *RF_I, double *RF_O, double *U_I, double *U_O,  
         int *lV, int *uV, int *iaS, int iD, double *F_I, double *F_O,
         int *nSizes ) 
  
  {
  int nValues = nSizes [ 3 - 1 ] * nSizes [ 2 - 1 ] * nSizes [ 1 - 1 ];
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x - 1 ) / block_Dim.x ) );
  
  iaS [ iD - 1 ] = -1;
  hipLaunchKernelGGL 
    ( ( ComputeFluxesInnerDeviceKernel ), grid_Dim, block_Dim, 0, 0,
      AP_I, AM_I, RF_I, RF_O, U_I, U_O, lV, uV, iaS, F_I, 
      nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] ); 
  
  DeviceSynchronize (  );
  
  iaS [ iD - 1 ] = +1;
  hipLaunchKernelGGL
    ( ( ComputeFluxesOuterDeviceKernel ), grid_Dim, block_Dim, 0, 0, 
      AP_O, AM_O, RF_I, RF_O, U_I, U_O, lV, uV, iaS, F_O, 
      nSizes [ 0 ], nSizes [ 1 ], nSizes [ 2 ] );
      
  DeviceSynchronize (  );
  }
                                                           
  
__global__ void ComputeUpdateDeviceKernel
                  ( double *dU, double *F_I, double *F_O,
                    double V, double A, double dT, int nValues )
  {
  if ( tiD <  nValues )
    dU [ tiD ] = dU [ tiD ] - dT * ( F_O [ tiD ] - F_I [ tiD ] ) * ( A / V );
  }


void ComputeUpdate_C
       ( double *dU, double *F_I, double *F_O,
         double V, double A, double dT, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x - 1 ) / block_Dim.x ) );
  
  hipLaunchKernelGGL
    ( ( ComputeUpdateDeviceKernel ), grid_Dim, block_Dim, 0, 0, 
      dU, F_I, F_O, V, A, dT, nValues );

  DeviceSynchronize ( );
  }


__global__ void AddUpdateDeviceKernel
                  ( double *O, double *U, double *C, int nValues )
  {
  if ( tiD < nValues )
    C [ tiD ] = O [ tiD ] + U [ tiD ];
  }
                

void AddUpdate_C ( double *O, double *U, double *C, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x - 1 ) / block_Dim.x ) );

  hipLaunchKernelGGL 
    ( ( AddUpdateDeviceKernel ), grid_Dim, block_Dim, 0, 0, 
      O, U, C, nValues );
    
  DeviceSynchronize ( );
  }


__global__ void CombineUpdatesDeviceKernel
                  ( double *C, double *O, double *U, int nValues )
  {
  if ( tiD < nValues )
    C [ tiD ] = 0.5 * ( O [ tiD ] + ( C [ tiD ] + U [ tiD ] ) );
  }
  
void CombineUpdates_C ( double *C, double *O, double *U, int nValues )
  {
  dim3 block_Dim ( BLOCK_DIM * BLOCK_DIM );
  dim3 grid_Dim  ( ceil ( ( nValues + block_Dim.x - 1 ) / block_Dim.x ) );
    
  hipLaunchKernelGGL 
    ( ( CombineUpdatesDeviceKernel ), grid_Dim, block_Dim, 0, 0, 
      C, O, U, nValues );
 
  DeviceSynchronize ( );
  }
