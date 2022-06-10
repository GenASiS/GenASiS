#include "Preprocessor"

submodule ( Slope_DFV_C_N__Form ) Slope_DFV_C_N__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeMomentumKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( M )
    
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then      
          S_S_1 ( iV )  =  - M ( iV )  *  N ( iV )  *  dPhi_1 ( iV )
          S_S_2 ( iV )  =  - M ( iV )  *  N ( iV )  *  dPhi_2 ( iV )
          S_S_3 ( iV )  =  - M ( iV )  *  N ( iV )  *  dPhi_3 ( iV )
        else
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then      
          S_S_1 ( iV )  =  - M ( iV )  *  N ( iV )  *  dPhi_1 ( iV )
          S_S_2 ( iV )  =  - M ( iV )  *  N ( iV )  *  dPhi_2 ( iV )
          S_S_3 ( iV )  =  - M ( iV )  *  N ( iV )  *  dPhi_3 ( iV )
        else
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeMomentumKernel


  module procedure ComputeEnergyKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( S_G )
    
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then      
          S_G ( iV )  =     V_1 ( iV )  *  S_S_1 ( iV )  &
                         +  V_2 ( iV )  *  S_S_2 ( iV )  &
                         +  V_3 ( iV )  *  S_S_3 ( iV )
        else
          S_G ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then      
          S_G ( iV )  =     V_1 ( iV )  *  S_S_1 ( iV )  &
                         +  V_2 ( iV )  *  S_S_2 ( iV )  &
                         +  V_3 ( iV )  *  S_S_3 ( iV )
        else
          S_G ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeEnergyKernel


end submodule Slope_DFV_C_N__Kernel
