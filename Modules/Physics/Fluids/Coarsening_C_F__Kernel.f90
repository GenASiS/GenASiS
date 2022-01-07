#include "Preprocessor"

submodule ( Coarsening_C_F__Form ) Coarsening_C_F__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iRZ, iPZ, &
      iR
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

      do iR  =  1, size ( nPZ )

        !$OMP OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET )
        do iPZ  =  1, nPZ ( iR )

          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + iPZ, :, iS_2 )  =  0.0_KDR
          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + iPZ, :, iS_3 )  =  0.0_KDR

          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_2 ) &
            =  0.0_KDR
          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_3 ) &
            =  0.0_KDR

        end do !-- iPZ
        !$OMP end OMP_TARGET_DIRECTIVE parallel do

      end do !-- iR

    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iRZ  =  1, nRZ
        FS_4D ( oC ( 1 ) + iRZ, :, :, iS_2 )  =  0.0_KDR
        FS_4D ( oC ( 1 ) + iRZ, :, :, iS_3 )  =  0.0_KDR
      end do !-- iRZ
      !$OMP end parallel do

      do iR  =  1, size ( nPZ )

        !$OMP parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST )
        do iPZ  =  1, nPZ ( iR )

          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + iPZ, :, iS_2 )  =  0.0_KDR
          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + iPZ, :, iS_3 )  =  0.0_KDR

          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_2 ) &
            =  0.0_KDR
          FS_4D ( oC ( 1 ) + iR, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_3 ) &
            =  0.0_KDR

        end do !-- iPZ
        !$OMP end parallel do

      end do !-- iR

    end if

  end procedure ComputeKernel


end submodule Coarsening_C_F__Kernel
