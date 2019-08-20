#include "Preprocessor"

submodule ( PressurelessFluid_Form ) PressurelessFluid_Kernel

  use Basics
  
  implicit none
    
contains


  module procedure ComputeConservedKernel

    integer ( KDI ) :: &
      iV
    
    if ( UseDevice ) then 
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, size ( D )
        D   ( iV ) = N ( iV )
        S_1 ( iV ) = N ( iV ) * V_1 ( iV )
        S_2 ( iV ) = N ( iV ) * V_2 ( iV )
        S_3 ( iV ) = N ( iV ) * V_3 ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else 
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, size ( D )
        D   ( iV ) = N ( iV )
        S_1 ( iV ) = N ( iV ) * V_1 ( iV )
        S_2 ( iV ) = N ( iV ) * V_2 ( iV )
        S_3 ( iV ) = N ( iV ) * V_3 ( iV )
      end do
      !$OMP end parallel do
    
    end if
    
  end procedure ComputeConservedKernel


  module procedure ComputePrimitiveKernel

    integer ( KDI ) :: &
      iV
      
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, size ( N )
        N ( iV ) = D ( iV )
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
    
    else
      
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, size ( N )
        N ( iV ) = D ( iV )
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
      !$OMP end parallel do
    
    end if
      
  end procedure ComputePrimitiveKernel


  module procedure ComputeEigenspeedsKernel

    integer ( KDI ) :: &
      iV
    
    if ( UseDevice ) then
    
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
      
    else 
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, size ( FEP_1 )
        FEP_1 ( iV ) = V_1 ( iV )
        FEP_2 ( iV ) = V_2 ( iV )
        FEP_3 ( iV ) = V_3 ( iV )
        FEM_1 ( iV ) = V_1 ( iV )
        FEM_2 ( iV ) = V_2 ( iV )
        FEM_3 ( iV ) = V_3 ( iV )
      end do
      !$OMP end parallel do
    
    end if
    
  end procedure ComputeEigenspeedsKernel


  module procedure ApplyBoundaryConditionsReflecting

    integer ( KDI ) :: &
      iV, jV, kV
      
    if ( UseDevice ) then
      
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
    
    else 
    
      !$OMP  parallel do collapse ( 3 ) &
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
      !$OMP end parallel do    
    
    end if
               
  end procedure ApplyBoundaryConditionsReflecting


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
    
    if ( UseDevice ) then
    
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
      
    else 
    
      !$OMP  parallel do collapse ( 3 ) &
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
      !$OMP end parallel do
      
    end if
      
  end procedure ComputeRiemannSolverInputKernel


end submodule PressurelessFluid_Kernel
