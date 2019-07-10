#include "Preprocessor"

submodule ( PolytropicFluid_Form ) PolytropicFluid_Kernel

  use iso_c_binding
  use Basics

  implicit none

contains


  module procedure ComputeConservedKernel

    G = E  +  0.5_KDR * N * ( V_1 **2  +  V_2 ** 2  +  V_3 ** 2 )

  end procedure ComputeConservedKernel


  module procedure ComputeConservedKernelDevice

    integer ( KDI ) :: &
      iV
      
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( G )
      G ( iV ) = E ( iV ) + 0.5_KDR * N ( iV ) &
                 * ( V_1 ( iV ) * V_1 ( iV ) &
                     + V_2 ( iV ) * V_2 ( iV ) &
                     + V_3 ( iV ) * V_3 ( iV ) )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    
  end procedure ComputeConservedKernelDevice


  module procedure ComputePrimitiveKernel

    real ( KDR ), dimension ( size ( N ) ) :: &
      KE
      
    !-- FIXME: 'associate' construct causes segfault with IBM XL
    !associate ( KE => 0.5_KDR * N * ( V_1 ** 2 + V_2 ** 2 + V_3 ** 2 ) )
    
    KE = 0.5_KDR * N * ( V_1 ** 2 + V_2 ** 2 + V_3 ** 2 )
    E  = G - KE
    where ( E < 0.0_KDR )
      E = 0.0_KDR
      G = KE
    end where

    !end associate !-- KE

  end procedure ComputePrimitiveKernel


  module procedure ComputePrimitiveKernelDevice

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      KE
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( KE )
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
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
  end procedure ComputePrimitiveKernelDevice


  module procedure ComputeAuxiliaryKernel

    P = E * ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end procedure ComputeAuxiliaryKernel


  module procedure ComputeAuxiliaryKernelDevice

    integer ( KDI ) :: &
      iV
      
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( P )
      P ( iV ) = E ( iV ) * ( Gamma ( iV ) - 1.0_KDR )
      if ( N ( iV ) > 0.0_KDR ) then
        K ( iV ) = P ( iV ) / ( N ( iV ) ** Gamma ( iV ) )
      else
        K ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
  end procedure ComputeAuxiliaryKernelDevice


  module procedure ComputeAuxiliaryFromPressureKernel

    E = P / ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end procedure ComputeAuxiliaryFromPressureKernel


  module procedure ComputeEigenspeedsKernel

    where ( N > 0.0_KDR .and. P > 0.0_KDR )
      CS = sqrt ( Gamma * P / N )
    elsewhere
      CS = 0.0_KDR
    end where 
    
    FEP_1 = V_1 + CS
    FEP_2 = V_2 + CS
    FEP_3 = V_3 + CS
    FEM_1 = V_1 - CS
    FEM_2 = V_2 - CS
    FEM_3 = V_3 - CS

  end procedure ComputeEigenspeedsKernel


  module procedure ComputeEigenspeedsKernelDevice

    integer ( KDI ) :: &
      iV
      
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( N )
      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
        CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / N ( iV ) )
      else
        CS ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( N )
      FEP_1 ( iV ) = V_1 ( iV ) + CS ( iV )
      FEP_2 ( iV ) = V_2 ( iV ) + CS ( iV )
      FEP_3 ( iV ) = V_3 ( iV ) + CS ( iV )
      FEM_1 ( iV ) = V_1 ( iV ) - CS ( iV )
      FEM_2 ( iV ) = V_2 ( iV ) - CS ( iV )
      FEM_3 ( iV ) = V_3 ( iV ) - CS ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
  end procedure ComputeEigenspeedsKernelDevice
  
  
  module procedure ApplyBoundaryConditionsReflecting

    E_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
          oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
          oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = E_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
              oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
              oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    Gamma_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
          oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
          oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = Gamma_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
              oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
              oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

  end procedure ApplyBoundaryConditionsReflecting


  module procedure ApplyBoundaryConditionsReflectingDevice

    integer ( KDI ) :: &
      iV, jV, kV
      
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE )
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
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
  end procedure ApplyBoundaryConditionsReflectingDevice


  module procedure ComputeRawFluxesKernel

    integer ( KDI ) :: &
      iV
      
    !F_D   = D   * V_Dim
    !F_S_1 = S_1 * V_Dim
    !F_S_2 = S_2 * V_Dim
    !F_S_3 = S_3 * V_Dim
    !F_S_Dim = F_S_Dim + P
    !F_G = ( G + P ) * V_Dim
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( F_D )
      F_D ( iV )     = D ( iV )   * V_Dim ( iV ) 
      F_S_1 ( iV )   = S_1 ( iV ) * V_Dim ( iV ) 
      F_S_2 ( iV )   = S_2 ( iV ) * V_Dim ( iV ) 
      F_S_3 ( iV )   = S_3 ( iV ) * V_Dim ( iV ) 
      F_S_Dim ( iV ) = F_S_Dim ( iV ) + P ( iV ) 
      F_G ( iV )     = ( G ( iV ) + P ( iV ) ) * V_Dim ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do

  end procedure ComputeRawFluxesKernel


end submodule PolytropicFluid_Kernel
