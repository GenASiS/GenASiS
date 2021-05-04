#include "Preprocessor"

submodule ( CurrentSet_C__Form ) CurrentSet_C__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeFluxesKernel

    integer ( KDI ) :: &
      iV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( F_D )
        F_D ( iV )  =  D ( iV )  *  V_Dim
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( F_D )
        F_D ( iV )  =  D ( iV )  *  V_Dim
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeFluxesKernel


  module procedure ComputeEigenspeedsKernel

    integer ( KDI ) :: &
      iV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, size ( EF_P )
        EF_P ( iV )  =  V_Dim
        EF_M ( iV )  =  V_Dim
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, size ( EF_P )
        EF_P ( iV )  =  V_Dim
        EF_M ( iV )  =  V_Dim
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeEigenspeedsKernel


end submodule CurrentSet_C__Kernel
