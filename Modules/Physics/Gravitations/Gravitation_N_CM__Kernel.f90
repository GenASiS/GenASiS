#include "Preprocessor"

submodule ( Gravitation_N_CM__Form ) Gravitation_N_CM__Kernel

  use Basics
  implicit none
  
contains


  module procedure SolveKernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
    logical ( KDL ) :: &
      UseDevice
    
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( Phi )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        Phi       ( iV )  =  - G * M  /  R ( iV )
        GradPhi_1 ( iV )  =    G * M  /  R ( iV ) ** 2
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        Phi       ( iV )  =  - G * M  /  R ( iV )
        GradPhi_1 ( iV )  =    G * M  /  R ( iV ) ** 2
      end do
      !$OMP end parallel do
    end if

  end procedure SolveKernel


end submodule Gravitation_N_CM__Kernel
