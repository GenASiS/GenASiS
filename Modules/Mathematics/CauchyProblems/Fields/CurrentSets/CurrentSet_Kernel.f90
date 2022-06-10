#include "Preprocessor"

submodule ( CurrentSet_Form ) CurrentSet_Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeEigenspeedsKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( EF_P )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        EF_P ( iV )  =  V_Dim ( iV )
        EF_M ( iV )  =  V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        EF_P ( iV )  =  V_Dim ( iV )
        EF_M ( iV )  =  V_Dim ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeEigenspeedsKernel


end submodule CurrentSet_Kernel
