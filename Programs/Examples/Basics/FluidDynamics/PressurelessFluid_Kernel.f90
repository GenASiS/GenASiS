#include "Preprocessor"

submodule ( PressurelessFluid_Form ) PressurelessFluid_Kernel

  use iso_c_binding
  use Basics
  
  implicit none
    
contains


  module procedure ComputeConservedKernel

    D   = N
    S_1 = N * V_1
    S_2 = N * V_2
    S_3 = N * V_3

  end procedure ComputeConservedKernel


  module procedure ComputeConservedKernelDevice

    integer ( KDI ) :: &
      iV
    
    call AssociateHost ( D_D, D )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( D )
      D   ( iV ) = N ( iV )
      S_1 ( iV ) = N ( iV ) * V_1 ( iV )
      S_2 ( iV ) = N ( iV ) * V_2 ( iV )
      S_3 ( iV ) = N ( iV ) * V_3 ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( D )

  end procedure ComputeConservedKernelDevice


  module procedure ComputePrimitiveKernel

    N = D

    where ( N > 0.0_KDR )
      V_1 = S_1 / N 
      V_2 = S_2 / N
      V_3 = S_3 / N
    elsewhere
      N   = 0.0_KDR
      V_1 = 0.0_KDR
      V_2 = 0.0_KDR
      V_3 = 0.0_KDR
      D   = 0.0_KDR
      S_1 = 0.0_KDR
      S_2 = 0.0_KDR
      S_3 = 0.0_KDR
    end where

  end procedure ComputePrimitiveKernel


  module procedure ComputePrimitiveKernelDevice

    integer ( KDI ) :: &
      iV

    call Copy ( D, D_D, D_N, N )
    
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_D, D )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( N )
      if ( N ( iV ) > 0.0_KDR ) then
        V_1 ( iV ) = S_1 ( iV ) / N ( iV )
        V_2 ( iV ) = S_2 ( iV ) / N ( iV )
        V_3 ( iV ) = S_3 ( iV ) / N ( iV )
      else
        N   ( iV )= 0.0_KDR
        V_1 ( iV ) = 0.0_KDR
        V_2 ( iV ) = 0.0_KDR
        V_3 ( iV ) = 0.0_KDR
        D   ( iV ) = 0.0_KDR
        S_1 ( iV ) = 0.0_KDR
        S_2 ( iV ) = 0.0_KDR
        S_3 ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( D )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    
  end procedure ComputePrimitiveKernelDevice

  
  module procedure ComputeEigenspeedsKernel

    FEP_1 = V_1
    FEP_2 = V_2
    FEP_3 = V_3
    FEM_1 = V_1
    FEM_2 = V_2
    FEM_3 = V_3

  end procedure ComputeEigenspeedsKernel


  module procedure ComputeEigenspeedsKernelDevice

    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_FEP_1, FEP_1 )
    call AssociateHost ( D_FEP_2, FEP_2 )
    call AssociateHost ( D_FEP_3, FEP_3 )
    call AssociateHost ( D_FEM_1, FEM_1 )
    call AssociateHost ( D_FEM_2, FEM_2 )
    call AssociateHost ( D_FEM_3, FEM_3 )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3,  V_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( FEP_1 )
      FEP_1 ( iV ) = V_1 ( iV )
      FEP_2 ( iV ) = V_2 ( iV )
      FEP_3 ( iV ) = V_3 ( iV )
      FEM_1 ( iV ) = V_1 ( iV )
      FEM_2 ( iV ) = V_2 ( iV )
      FEM_3 ( iV ) = V_3 ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( FEM_3 )
    call DisassociateHost ( FEM_2 )
    call DisassociateHost ( FEM_1 )
    call DisassociateHost ( FEP_3 )
    call DisassociateHost ( FEP_2 )
    call DisassociateHost ( FEP_1 )

  end procedure ComputeEigenspeedsKernelDevice


  module procedure ApplyBoundaryConditionsReflecting

    N_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
          oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
          oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = N_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
              oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
              oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    VI_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
           oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
           oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = - VI_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
                 oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
                 oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    VJ_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
           oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
           oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = VJ_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
               oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
               oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    VK_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
           oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
           oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = VK_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
               oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
               oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

  end procedure ApplyBoundaryConditionsReflecting


  module procedure ApplyBoundaryConditionsReflectingDevice

    integer ( KDI ) :: &
      iV, jV, kV
      
    call AssociateHost ( D_N_E, N_E )
    call AssociateHost ( D_VI_E, VI_E )
    call AssociateHost ( D_VJ_E, VJ_E )
    call AssociateHost ( D_VK_E, VK_E )
    call AssociateHost ( D_N_I, N_I )
    call AssociateHost ( D_VI_I, VI_I )
    call AssociateHost ( D_VJ_I, VJ_I )
    call AssociateHost ( D_VK_I, VK_I )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE )
    do kV = 1, nB ( 3 )
      do jV = 1, nB ( 2 )
        do iV = 1, nB ( 1 )

          N_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = N_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

          VI_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = - VI_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

          VJ_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = VJ_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

          VK_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = VK_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
               
    call DisassociateHost ( VK_I )
    call DisassociateHost ( VJ_I )
    call DisassociateHost ( VI_I )
    call DisassociateHost ( N_I )
    call DisassociateHost ( VK_E )
    call DisassociateHost ( VJ_E )
    call DisassociateHost ( VI_E )
    call DisassociateHost ( N_E )
    
  end procedure ApplyBoundaryConditionsReflectingDevice


  module procedure ComputeRawFluxesKernel

    F_D   = D   * V_Dim
    F_S_1 = S_1 * V_Dim
    F_S_2 = S_2 * V_Dim
    F_S_3 = S_3 * V_Dim

  end procedure ComputeRawFluxesKernel


  module procedure ComputeRiemannSolverInputKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS_P, iaS_M, &
      iaVS_P, iaVS_M, &
      lV, uV
      
    !AP_I = max ( 0.0_KDR, + cshift ( LP_O, shift = -1, dim = iD ), + LP_I )
    !AP_O = max ( 0.0_KDR, + LP_O, + cshift ( LP_I, shift = +1, dim = iD ) )
    !AM_I = max ( 0.0_KDR, - cshift ( LM_O, shift = -1, dim = iD ), - LM_I )
    !AM_O = max ( 0.0_KDR, - LM_O, - cshift ( LM_I, shift = +1, dim = iD ) )
    
    call AssociateHost ( D_AP_I, AP_I )
    call AssociateHost ( D_AP_O, AP_O )
    call AssociateHost ( D_AM_I, AM_I )
    call AssociateHost ( D_AM_O, AM_O )
    call AssociateHost ( D_LP_I, LP_I )
    call AssociateHost ( D_LP_O, LP_O )
    call AssociateHost ( D_LM_I, LM_I )
    call AssociateHost ( D_LM_O, LM_O )

    lV = 1
    where ( shape ( LP_O ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV
    
    uV = 1
    where ( shape ( LP_O ) > 1 )
      uV = shape ( LP_O ) - oV
    end where
    uV ( iD ) = size ( LP_O, dim = iD ) - 1
    
    iaS_M = 0
    iaS_M ( iD ) = - 1
    iaS_P = 0
    iaS_P ( iD ) = + 1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS_M, iaVS_P )
    do kV = lV ( 3 ), uV ( 3 )    
      do jV = lV ( 2 ), uV ( 2 )  
        do iV = lV ( 1 ), uV ( 1 )
        
          iaVS_M = [ iV, jV, kV ] + iaS_M
          iaVS_P = [ iV, jV, kV ] + iaS_P
          
          AP_I ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    + LP_O ( iaVS_M ( 1 ), iaVS_M ( 2 ), iaVS_M ( 3 ) ), &
                    + LP_I ( iV, jV, kV ) )
          AP_O ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    + LP_O ( iV, jV, kV ), &
                    + LP_I ( iaVS_P ( 1 ), iaVS_P ( 2 ), iaVS_P ( 3 ) ) )
          AM_I ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    - LM_O ( iaVS_M ( 1 ), iaVS_M ( 2 ), iaVS_M ( 3 ) ), &
                    - LM_I ( iV, jV, kV ) )
          AM_O ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    - LM_O ( iV, jV, kV ), &
                    - LM_I ( iaVS_P ( 1 ), iaVS_P ( 2 ), iaVS_P ( 3 ) ) )
        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( LM_O )
    call DisassociateHost ( LM_I )
    call DisassociateHost ( LP_O )
    call DisassociateHost ( LP_I )
    call DisassociateHost ( AM_O )
    call DisassociateHost ( AM_I )
    call DisassociateHost ( AP_O )
    call DisassociateHost ( AP_I )
    
  end procedure ComputeRiemannSolverInputKernel


end submodule PressurelessFluid_Kernel