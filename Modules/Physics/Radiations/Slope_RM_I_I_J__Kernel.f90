#include "Preprocessor"

submodule ( Slope_RM_I_I_J__Form ) Slope_RM_I_I_J__Kernel
  
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
        else
          S_E   ( iV )  =  0.0_KDR
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
          S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) ) &
                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
        else
          S_E   ( iV )  =  0.0_KDR
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


end submodule Slope_RM_I_I_J__Kernel
