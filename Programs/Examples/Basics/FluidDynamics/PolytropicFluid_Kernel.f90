#include "Preprocessor"

submodule ( PolytropicFluid_Form ) PolytropicFluid_Kernel

  use iso_c_binding
  use Basics
  
  implicit none
  
  include 'PolytropicFluid_Interface.f90'

contains


  module procedure ComputeConservedKernel

    integer ( KDI ) :: &
      iV, &
      nV
      
    nV = size ( G )

    if ( UseDevice ) then
          
      if ( UseDirectDevice ) then
        !$OMP target data use_device_ptr &
        !$OMP  ( G, E, N, V_1, V_2, V_3 )    
        call ComputeConservedPolytropic_C &
               ( c_loc ( G ), c_loc ( E ), c_loc ( N ), &
                 c_loc ( V_1 ), c_loc ( V_2 ), c_loc ( V_3 ), nV )
        !$OMP end target data
        return
      end if

      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( G )
        G ( iV ) = E ( iV ) + 0.5_KDR * N ( iV ) &
                   * ( V_1 ( iV ) * V_1 ( iV ) &
                       + V_2 ( iV ) * V_2 ( iV ) &
                       + V_3 ( iV ) * V_3 ( iV ) )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
  
    else      

      !$OMP  parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
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
      iV, &
      nV 
    real ( KDR ) :: &
      KE
    
    nV = size ( E )

    if ( UseDevice ) then
    
      if ( UseDirectDevice ) then
        !$OMP target data use_device_ptr &
        !$OMP   ( E, G, N, V_1, V_2, V_3 )    
        call ComputePrimitivePolytropic_C &
               ( c_loc ( E ), c_loc ( G ), c_loc ( N ), &
                 c_loc ( V_1 ), c_loc ( V_2 ), c_loc ( V_3 ), nV )
        !$OMP end target data
        return
      end if
    
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( KE )
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
    
      !$OMP  parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( KE )
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
      iV, &
      nV
      
    nV = size ( P )
      
    if ( UseDevice ) then
    
      if ( UseDirectDevice ) then
        !$OMP target data use_device_ptr &
        !$OMP   ( P, K, N, E, Gamma )    
        call ComputeAuxiliaryPolytropic_C &
             ( c_loc ( P ), c_loc ( K ), c_loc ( N ), &
               c_loc ( E ), c_loc ( Gamma ), nV )
        !$OMP end target data
        return
      end if
      
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
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

      !$OMP  parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
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
      iV, &
      nV
      
    nV = size ( N )
    
    if ( UseDevice ) then
    
      if ( UseDirectDevice ) then
        !$OMP target data use_device_ptr &
        !$OMP   ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
        !$OMP     CS, N, V_1, V_2, V_3, P, Gamma ) 
        call ComputeEigenspeedsPolytropic_C &
               ( c_loc ( FEP_1 ), c_loc ( FEP_2 ), c_loc ( FEP_3 ), &
                 c_loc ( FEM_1 ), c_loc ( FEM_2 ), c_loc ( FEM_3 ), &
                 c_loc ( CS ), c_loc ( N ), &
                 c_loc ( V_1 ), c_loc ( V_2 ), c_loc ( V_3 ), &
                 c_loc ( P ), c_loc ( Gamma ), nV )
        !$OMP end target data
        return
      end if
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( N )
        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / N ( iV ) )
        else
          CS ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
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
    
      !$OMP  parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( N )
        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / N ( iV ) )
        else
          CS ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end parallel do simd
      
      !$OMP  parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
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
    
    integer ( KDI ), dimension ( 3 ) :: &
      nSizes
      
    nSizes = shape ( E_E )
      
    if ( UseDevice ) then
    
      if ( UseDirectDevice ) then
        !$OMP target data use_device_ptr ( E_E, Gamma_E, E_I, Gamma_I )
        call ApplyBoundaryConditionsReflectingPolytropic_C &
               ( c_loc ( E_E ), c_loc ( Gamma_E ), & 
                 c_loc ( E_I ), c_loc ( Gamma_I ), &
                 c_loc ( nB ), c_loc ( oBE ), c_loc ( oBI ), nSizes )
        !$OMP end target data
        return
      end if

      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
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
    
      !$OMP  parallel do simd collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
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

    integer ( KDI ) :: &
      iV, &
      nV
      
    !F_D   = D   * V_Dim
    !F_S_1 = S_1 * V_Dim
    !F_S_2 = S_2 * V_Dim
    !F_S_3 = S_3 * V_Dim
    !F_S_Dim = F_S_Dim + P
    !F_G = ( G + P ) * V_Dim
    
    nV = size ( F_D )
    
    if ( UseDevice ) then
    
      if ( UseDirectDevice ) then
        !$OMP target data use_device_ptr &
        !$OMP   ( F_D, F_S_1, F_S_2, F_S_3, &
        !$OMP     F_S_Dim, F_G, D, S_1, S_2, S_3, G, P, V_Dim )
        call ComputeRawFluxesPolytropic_C &
               ( c_loc ( F_D ), &
                 c_loc ( F_S_1 ), c_loc ( F_S_2 ), c_loc ( F_S_3 ), &
                 c_loc ( F_S_Dim ), c_loc ( F_G ), c_loc ( D ), &
                 c_loc ( S_1 ), c_loc ( S_2 ), c_loc ( S_3 ), &
                 c_loc ( G ), c_loc ( P ), c_loc ( V_Dim ), nV )
        !$OMP  end target data
        return
      end if
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( F_D )
        F_D ( iV )     = D ( iV )   * V_Dim ( iV ) 
        F_S_1 ( iV )   = S_1 ( iV ) * V_Dim ( iV ) 
        F_S_2 ( iV )   = S_2 ( iV ) * V_Dim ( iV ) 
        F_S_3 ( iV )   = S_3 ( iV ) * V_Dim ( iV ) 
        F_S_Dim ( iV ) = F_S_Dim ( iV ) + P ( iV ) 
        F_G ( iV )     = ( G ( iV ) + P ( iV ) ) * V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
    else
      
      !$OMP  parallel do simd &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( F_D )
        F_D ( iV )     = D ( iV )   * V_Dim ( iV ) 
        F_S_1 ( iV )   = S_1 ( iV ) * V_Dim ( iV ) 
        F_S_2 ( iV )   = S_2 ( iV ) * V_Dim ( iV ) 
        F_S_3 ( iV )   = S_3 ( iV ) * V_Dim ( iV ) 
        F_S_Dim ( iV ) = F_S_Dim ( iV ) + P ( iV ) 
        F_G ( iV )     = ( G ( iV ) + P ( iV ) ) * V_Dim ( iV )
      end do
      !$OMP end parallel do simd
      
    end if
      
  end procedure ComputeRawFluxesKernel


end submodule PolytropicFluid_Kernel
