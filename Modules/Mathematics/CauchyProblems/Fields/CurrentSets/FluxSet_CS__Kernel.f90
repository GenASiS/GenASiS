#include "Preprocessor"

submodule ( FluxSet_CS__Form ) FluxSet_CS__Kernel
  
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
      
    nV  =  size ( F_D )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


end submodule FluxSet_CS__Kernel
