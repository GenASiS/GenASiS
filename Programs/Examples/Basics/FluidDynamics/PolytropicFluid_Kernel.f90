#include "Preprocessor"

submodule ( PolytropicFluid_Form ) PolytropicFluid_Kernel

  use Basics
  
  implicit none

contains


  module procedure ComputeConservedKernel

    integer ( KDI ) :: &
      iV
      
    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( G )
        G ( iV ) = E ( iV ) + 0.5_KDR * N ( iV ) &
                   * ( V_1 ( iV ) * V_1 ( iV ) &
                       + V_2 ( iV ) * V_2 ( iV ) &
                       + V_3 ( iV ) * V_3 ( iV ) )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
  
    else      

      !$OMP parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( G )
        G ( iV ) = E ( iV ) + 0.5_KDR * N ( iV ) &
                   * ( V_1 ( iV ) * V_1 ( iV ) &
                       + V_2 ( iV ) * V_2 ( iV ) &
                       + V_3 ( iV ) * V_3 ( iV ) )
      end do
      !$OMP end parallel do simd
      
    end if
    
    
  end procedure ComputeConservedKernel


  module procedure ComputePrimitiveKernel

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      KE
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( KE )
      do iV = 1, size ( E )
      
        KE = 0.5_KDR * N ( iV ) &
               * ( V_1 ( iV ) * V_1 ( iV )  +  V_2 ( iV ) * V_2 ( iV ) &
                   + V_3 ( iV ) * V_3 ( iV ) )
        E ( iV )  = G ( iV ) - KE

        if ( E ( iV ) < 0.0_KDR ) then
          E ( iV ) = 0.0_KDR
          G ( iV ) = KE
        end if

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
      
    else 
    
      !$OMP parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( KE )
      do iV = 1, size ( E )
      
        KE = 0.5_KDR * N ( iV ) &
               * ( V_1 ( iV ) * V_1 ( iV )  +  V_2 ( iV ) * V_2 ( iV ) &
                   + V_3 ( iV ) * V_3 ( iV ) )
        E ( iV )  = G ( iV ) - KE

        if ( E ( iV ) < 0.0_KDR ) then
          E ( iV ) = 0.0_KDR
          G ( iV ) = KE
        end if

      end do
      !$OMP end parallel do simd
      
    end if
      
  end procedure ComputePrimitiveKernel


  module procedure ComputeAuxiliaryKernel

    integer ( KDI ) :: &
      iV
      
    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( P )
        P ( iV ) = E ( iV ) * ( Gamma ( iV ) - 1.0_KDR )
        if ( N ( iV ) > 0.0_KDR ) then
          K ( iV ) = P ( iV ) / ( N ( iV ) ** Gamma ( iV ) )
        else
          K ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
    else      

      !$OMP parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( P )
        P ( iV ) = E ( iV ) * ( Gamma ( iV ) - 1.0_KDR )
        if ( N ( iV ) > 0.0_KDR ) then
          K ( iV ) = P ( iV ) / ( N ( iV ) ** Gamma ( iV ) )
        else
          K ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end parallel do simd
    
    end if
      
  end procedure ComputeAuxiliaryKernel


  module procedure ComputeAuxiliaryFromPressureKernel

    E = P / ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end procedure ComputeAuxiliaryFromPressureKernel


  module procedure ComputeEigenspeedsKernel

    integer ( KDI ) :: &
      iV
    
    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( N )
        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / N ( iV ) )
        else
          CS ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( N )
        FEP_1 ( iV ) = V_1 ( iV ) + CS ( iV )
        FEP_2 ( iV ) = V_2 ( iV ) + CS ( iV )
        FEP_3 ( iV ) = V_3 ( iV ) + CS ( iV )
        FEM_1 ( iV ) = V_1 ( iV ) - CS ( iV )
        FEM_2 ( iV ) = V_2 ( iV ) - CS ( iV )
        FEM_3 ( iV ) = V_3 ( iV ) - CS ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
      
    else
    
      !$OMP parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( N )
        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / N ( iV ) )
        else
          CS ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end parallel do simd
      
      !$OMP parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( N )
        FEP_1 ( iV ) = V_1 ( iV ) + CS ( iV )
        FEP_2 ( iV ) = V_2 ( iV ) + CS ( iV )
        FEP_3 ( iV ) = V_3 ( iV ) + CS ( iV )
        FEM_1 ( iV ) = V_1 ( iV ) - CS ( iV )
        FEM_2 ( iV ) = V_2 ( iV ) - CS ( iV )
        FEM_3 ( iV ) = V_3 ( iV ) - CS ( iV )
      end do
      !$OMP end parallel do simd
    
    end if
    
  end procedure ComputeEigenspeedsKernel
  
  
  module procedure ApplyBoundaryConditionsReflecting

    integer ( KDI ) :: &
      iV, jV, kV
      
    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
          
            E_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
              = E_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )
            
            Gamma_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
              = Gamma_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )
              
          end do
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
      
    else
    
      !$OMP parallel do simd collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
          
            E_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
              = E_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )
            
            Gamma_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
              = Gamma_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )
              
          end do
        end do
      end do
      !$OMP end parallel do simd
    
    end if
    
  end procedure ApplyBoundaryConditionsReflecting


  module procedure ComputeRawFluxesKernel

    integer :: &
      Delta_1, Delta_2, Delta_3
    integer ( KDI ) :: &
      iV
      
    Delta_1  =  int (    real ( 1 + iDim  -  abs ( 1 - iDim ) )  &
                      /  real ( 1 + iDim  +  abs ( 1 - iDim ) ) )

    Delta_2  =  int (    real ( 2 + iDim  -  abs ( 2 - iDim ) )  &
                      /  real ( 2 + iDim  +  abs ( 2 - iDim ) ) )

    Delta_3  =  int (    real ( 3 + iDim  -  abs ( 3 - iDim ) )  &
                      /  real ( 3 + iDim  +  abs ( 3 - iDim ) ) )
      
    !F_D   = D   * V_Dim
    !F_S_1 = S_1 * V_Dim
    !F_S_2 = S_2 * V_Dim
    !F_S_3 = S_3 * V_Dim
    !F_S_Dim = F_S_Dim + P
    !F_G = ( G + P ) * V_Dim
    
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( P )
        F_D ( iV )     = D ( iV )   * V_Dim ( iV ) 
        F_S_1 ( iV )   = S_1 ( iV ) * V_Dim ( iV )  +  Delta_1 * P ( iV ) 
        F_S_2 ( iV )   = S_2 ( iV ) * V_Dim ( iV )  +  Delta_2 * P ( iV ) 
        F_S_3 ( iV )   = S_3 ( iV ) * V_Dim ( iV )  +  Delta_3 * P ( iV ) 
        F_G ( iV )     = ( G ( iV ) + P ( iV ) ) * V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
    else
      !$OMP parallel do simd &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( P )
        F_D ( iV )     = D ( iV )   * V_Dim ( iV ) 
        F_S_1 ( iV )   = S_1 ( iV ) * V_Dim ( iV )  +  Delta_1 * P ( iV ) 
        F_S_2 ( iV )   = S_2 ( iV ) * V_Dim ( iV )  +  Delta_2 * P ( iV ) 
        F_S_3 ( iV )   = S_3 ( iV ) * V_Dim ( iV )  +  Delta_3 * P ( iV ) 
        F_G ( iV )     = ( G ( iV ) + P ( iV ) ) * V_Dim ( iV )
      end do
      !$OMP end parallel do simd
    end if
      
  end procedure ComputeRawFluxesKernel


end submodule PolytropicFluid_Kernel
