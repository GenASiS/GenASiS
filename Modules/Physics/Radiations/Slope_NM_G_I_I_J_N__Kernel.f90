#include "Preprocessor"

submodule ( Slope_NM_G_I_I_J_N__Form ) Slope_NM_G_I_I_J_N__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( S_E )
    
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then      

          S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) ) &
                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
          S_D   ( iV )  =  ( Xi_N ( iV )  -  Chi_N ( iV )  *  N ( iV ) ) &
                           /  ( 1.0_KDR  +  Chi_N ( iV ) * dT )

          ! !-- Limiter
          ! S_E ( iV )  &
          !   =  sign ( min ( abs ( S_E ( iV ) ), &
          !                   max ( 1.e-2 * abs ( J ( iV ) ) / dT, &
          !                         1.e-4 * abs ( E_F ( iV ) ) / dT ) ), &
          !             S_E ( iV ) )
          ! S_D ( iV )  &
          !   =  sign ( min ( abs ( S_D ( iV ) ), &
          !                   max ( 1.e-2 * abs ( N ( iV ) ) / dT, &
          !                         1.e-4 * abs ( Y_F ( iV ) &
          !                                       * N_F ( iV ) ) / dT ) ), &
          !             S_D ( iV ) )

        else
          S_E   ( iV )  =  0.0_KDR
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
          S_D   ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then      

          S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) ) &
                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
          S_D   ( iV )  =  ( Xi_N ( iV )  -  Chi_N ( iV )  *  N ( iV ) ) &
                           /  ( 1.0_KDR  +  Chi_N ( iV ) * dT )

          ! !-- Limiter
          ! S_E ( iV )  &
          !   =  sign ( min ( abs ( S_E ( iV ) ), &
          !                   max ( 1.e-2 * abs ( J ( iV ) ) / dT, &
          !                         1.e-4 * abs ( E_F ( iV ) ) / dT ) ), &
          !             S_E ( iV ) )
          ! S_D ( iV )  &
          !   =  sign ( min ( abs ( S_D ( iV ) ), &
          !                   max ( 1.e-2 * abs ( N ( iV ) ) / dT, &
          !                         1.e-4 * abs ( Y_F ( iV ) &
          !                                       * N_F ( iV ) ) / dT ) ), &
          !             S_D ( iV ) )

        else
          S_E   ( iV )  =  0.0_KDR
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
          S_D   ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


end submodule Slope_NM_G_I_I_J_N__Kernel
