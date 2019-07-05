#include "Preprocessor"

module PolytropicFluid_Kernel

  use iso_c_binding
  use Basics

  implicit none
  private
  
  public :: &
    ComputeConservedKernel, &
    ComputeConservedKernelDevice, &
    ComputePrimitiveKernel, &
    ComputePrimitiveKernelDevice, &
    ComputeAuxiliaryKernel, &
    ComputeAuxiliaryKernelDevice, &
    ComputeAuxiliaryFromPressureKernel, &
    ComputeEigenspeedsKernel, &
    ComputeEigenspeedsKernelDevice, &
    ApplyBoundaryConditionsReflecting, &
    ApplyBoundaryConditionsReflectingDevice, &
    ComputeRawFluxesKernel
    
contains


  subroutine ComputeConservedKernel ( G, E, N, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      E, &
      N, &
      V_1, V_2, V_3

    G = E  +  0.5_KDR * N * ( V_1 **2  +  V_2 ** 2  +  V_3 ** 2 )

  end subroutine ComputeConservedKernel


  subroutine ComputeConservedKernelDevice &
               ( G, E, N, V_1, V_2, V_3, &
                 D_G, D_E, D_N, D_V_1, D_V_2, D_V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      E, &
      N, &
      V_1, V_2, V_3
    type ( c_ptr ), intent ( in ) :: &
      D_G, &
      D_E, &
      D_N, &
      D_V_1, D_V_2, D_V_3
      
    integer ( KDI ) :: &
      iV
    
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( G )
      G ( iV ) = E ( iV ) + 0.5_KDR * N ( iV ) &
                 * ( V_1 ( iV ) * V_1 ( iV ) &
                     + V_2 ( iV ) * V_2 ( iV ) &
                     + V_3 ( iV ) * V_3 ( iV ) )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    
  end subroutine ComputeConservedKernelDevice


  subroutine ComputePrimitiveKernel ( E, G, N, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3
      
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

  end subroutine ComputePrimitiveKernel


  subroutine ComputePrimitiveKernelDevice &
               ( E, G, N, V_1, V_2, V_3, &
                 D_E, D_G, D_N, D_V_1, D_V_2, D_V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3
    type ( c_ptr ), intent ( in ) :: &
      D_E, &
      D_G, &
      D_N, &
      D_V_1, D_V_2, D_V_3
      
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      KE
    
    
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    
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
    
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( G )
    call DisassociateHost ( E )


  end subroutine ComputePrimitiveKernelDevice


  subroutine ComputeAuxiliaryKernel ( P, K, N, E, Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      E, &
      Gamma

    P = E * ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end subroutine ComputeAuxiliaryKernel


  subroutine ComputeAuxiliaryKernelDevice &
               ( P, K, N, E, Gamma, &
                 D_P, D_K, D_N, D_E, D_Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      E, &
      Gamma
    type ( c_ptr ), intent ( in ) :: &
      D_P, &
      D_K, &
      D_N, &
      D_E, &
      D_Gamma
    
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_K, K )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_Gamma, Gamma )
    
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
      
    call DisassociateHost ( Gamma )
    call DisassociateHost ( E )
    call DisassociateHost ( N )
    call DisassociateHost ( K )
    call DisassociateHost ( P )

  end subroutine ComputeAuxiliaryKernelDevice


  subroutine ComputeAuxiliaryFromPressureKernel ( E, K, N, P, Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      P, &
      Gamma

    E = P / ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end subroutine ComputeAuxiliaryFromPressureKernel


  subroutine ComputeEigenspeedsKernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, N, &
                 V_1, V_2, V_3, P, Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3, &
      P, &
      Gamma

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

  end subroutine ComputeEigenspeedsKernel


  subroutine ComputeEigenspeedsKernelDevice &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, N, &
                 V_1, V_2, V_3, P, Gamma, &
                 D_FEP_1, D_FEP_2, D_FEP_3, D_FEM_1, D_FEM_2, D_FEM_3, &
                 D_CS, D_N, D_V_1, D_V_2, D_V_3, D_P, D_Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3, &
      P, &
      Gamma
    type ( c_ptr ), intent ( in ) :: &
      D_FEP_1, D_FEP_2, D_FEP_3, &
      D_FEM_1, D_FEM_2, D_FEM_3, &
      D_CS, &
      D_N, &
      D_V_1, D_V_2, D_V_3, &
      D_P, &
      D_Gamma
      
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_FEP_1, FEP_1 )
    call AssociateHost ( D_FEP_2, FEP_2 )
    call AssociateHost ( D_FEP_3, FEP_3 )
    call AssociateHost ( D_FEM_1, FEM_1 )
    call AssociateHost ( D_FEM_2, FEM_2 )
    call AssociateHost ( D_FEM_3, FEM_3 )
    call AssociateHost ( D_CS, CS )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_Gamma, Gamma )
    
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
    
    call DisassociateHost ( Gamma )
    call DisassociateHost ( P )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( CS )
    call DisassociateHost ( FEM_3 )
    call DisassociateHost ( FEM_2 )
    call DisassociateHost ( FEM_1 )
    call DisassociateHost ( FEP_3 )
    call DisassociateHost ( FEP_2 )
    call DisassociateHost ( FEP_1 )

  end subroutine ComputeEigenspeedsKernelDevice
  
  
  subroutine ApplyBoundaryConditionsReflecting &
               ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      E_E, &
      Gamma_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      E_I, &
      Gamma_I
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB,  & 
      oBE, &
      oBI

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

  end subroutine ApplyBoundaryConditionsReflecting


  subroutine ApplyBoundaryConditionsReflectingDevice &
               ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI, &
                 D_E_E, D_Gamma_E, D_E_I, D_Gamma_I )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      E_E, &
      Gamma_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      E_I, &
      Gamma_I
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB,  & 
      oBE, &
      oBI
    type ( c_ptr ), intent ( in ) :: &
      D_E_E, &
      D_Gamma_E, &
      D_E_I, &
      D_Gamma_I
      
    integer ( KDI ) :: &
      iV, jV, kV
      
    call AssociateHost ( D_E_E, E_E )
    call AssociateHost ( D_Gamma_E, Gamma_E )
    call AssociateHost ( D_E_I, E_I )
    call AssociateHost ( D_Gamma_I, Gamma_I )
    
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
    
    call DisassociateHost ( Gamma_I )
    call DisassociateHost ( E_I )
    call DisassociateHost ( Gamma_E )
    call DisassociateHost ( E_E )

  end subroutine ApplyBoundaryConditionsReflectingDevice


  subroutine ComputeRawFluxesKernel &
               ( F_D, F_S_1, F_S_2, F_S_3, F_S_Dim, F_G, D, S_1, S_2, S_3, &
                 G, P, V_Dim, D_F_D, D_F_S_1, D_F_S_2, D_F_S_3, D_F_S_Dim, &
                 D_F_G, D_D, D_S_1, D_S_2, D_S_3, D_G, D_P, D_V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_D, &
      F_S_1, F_S_2, F_S_3, &
      F_S_Dim, &
      F_G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      D, &
      S_1, S_2, S_3, &
      G, &
      P, &
      V_Dim
    type ( c_ptr ), intent ( in ) :: &
      D_F_D, &
      D_F_S_1, D_F_S_2, D_F_S_3, &
      D_F_S_Dim, &
      D_F_G, &
      D_D, &
      D_S_1, D_S_2, D_S_3, &
      D_G, &
      D_P, &
      D_V_Dim
      
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_F_D, F_D )
    call AssociateHost ( D_F_S_1, F_S_1 )
    call AssociateHost ( D_F_S_2, F_S_2 )
    call AssociateHost ( D_F_S_3, F_S_3 )
    call AssociateHost ( D_F_S_Dim, F_S_Dim )
    call AssociateHost ( D_F_G, F_G )
    call AssociateHost ( D_D, D )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_V_Dim, V_Dim )
    
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

    call DisassociateHost ( V_Dim )
    call DisassociateHost ( P )
    call DisassociateHost ( G )
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( D )
    call DisassociateHost ( F_G )
    call DisassociateHost ( F_S_Dim )
    call DisassociateHost ( F_S_3 )
    call DisassociateHost ( F_S_2 )
    call DisassociateHost ( F_S_1 )
    call DisassociateHost ( F_D )
    
  end subroutine ComputeRawFluxesKernel


end module PolytropicFluid_Kernel
