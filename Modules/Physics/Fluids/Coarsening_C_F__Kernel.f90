#include "Preprocessor"

submodule ( Coarsening_C_F__Form ) Coarsening_C_F__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iRZ
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iRZ  =  1, nRZ

        FS_4D ( oC ( 1 ) + iRZ, :, :, iS_2 )  =  0.0_KDR
        FS_4D ( oC ( 1 ) + iRZ, :, :, iS_3 )  =  0.0_KDR

      end do !-- iRZ
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iRZ  =  1, nRZ

        FS_4D ( oC ( 1 ) + iRZ, :, :, iS_2 )  =  0.0_KDR
        FS_4D ( oC ( 1 ) + iRZ, :, :, iS_3 )  =  0.0_KDR

      end do !-- iRZ
      !$OMP end parallel do

    end if

  end procedure ComputeKernel


end submodule Coarsening_C_F__Kernel
