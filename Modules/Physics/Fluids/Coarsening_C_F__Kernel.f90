#include "Preprocessor"

submodule ( Coarsening_C_F__Form ) Coarsening_C_F__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iRZ, iPZ
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    if ( UseDevice ) then

      !-- First radial shell only
      if ( iaB ( 1 )  ==  1 ) then
        !$OMP OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET )
        do iRZ  =  1, nRZ
          FS_4D ( oC ( 1 ) + iRZ, :, :, iS_2 )  =  0.0_KDR
          FS_4D ( oC ( 1 ) + iRZ, :, :, iS_3 )  =  0.0_KDR
        end do !-- iRZ
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      end if

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iPZ  =  1, nPZ

        FS_4D ( :, oC ( 2 ) + iPZ, :, iS_2 )  =  0.0_KDR
        FS_4D ( :, oC ( 2 ) + iPZ, :, iS_3 )  =  0.0_KDR

        FS_4D ( :, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_2 )  =  0.0_KDR
        FS_4D ( :, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_3 )  =  0.0_KDR

      end do !-- iPZ
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !-- First radial shell only
      if ( iaB ( 1 )  ==  1 ) then
        !$OMP OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET )
        do iRZ  =  1, nRZ
          FS_4D ( oC ( 1 ) + iRZ, :, :, iS_2 )  =  0.0_KDR
          FS_4D ( oC ( 1 ) + iRZ, :, :, iS_3 )  =  0.0_KDR
        end do !-- iRZ
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      end if

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iPZ  =  1, nPZ

        FS_4D ( :, oC ( 2 ) + iPZ, :, iS_2 )  =  0.0_KDR
        FS_4D ( :, oC ( 2 ) + iPZ, :, iS_3 )  =  0.0_KDR

        FS_4D ( :, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_2 )  =  0.0_KDR
        FS_4D ( :, oC ( 2 ) + nC ( 2 ) - ( iPZ - 1 ), :, iS_3 )  =  0.0_KDR

      end do !-- iPZ
      !$OMP end parallel do

    end if

  end procedure ComputeKernel


  module procedure ComputeMoreKernel 

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( FS, dim = 1 )
    
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        if (      ( BP ( iV )  >=  0.99  .and. BP ( iV )  <  0.99 * nBT ) &
             .or. ( BA ( iV )  >=  0.99  .and. BA ( iV )  <  0.99 * nBT ) ) &
        then
          FS ( iV, iS_2 )  =  0.0_KDR
          FS ( iV, iS_3 )  =  0.0_KDR
        end if
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        if (      ( BP ( iV )  >=  0.99  .and. BP ( iV )  <  0.99 * nBT ) &
             .or. ( BA ( iV )  >=  0.99  .and. BA ( iV )  <  0.99 * nBT ) ) &
        then
          FS ( iV, iS_2 )  =  0.0_KDR
          FS ( iV, iS_3 )  =  0.0_KDR
        end if
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeMoreKernel


end submodule Coarsening_C_F__Kernel
