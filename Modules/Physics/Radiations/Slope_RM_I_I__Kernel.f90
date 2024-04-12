#include "Preprocessor"

submodule ( Slope_RM_I_I__Form ) Slope_RM_I_I__Kernel
  
  use Basics
  
  implicit none
  
contains


!   module procedure ComputeKernel

!     integer ( KDI ) :: &
!       iV, &
!       nV
!     logical ( KDL ) :: &
!       UseDevice

!     UseDevice = .false.
!     if ( present ( UseDeviceOption ) ) &
!       UseDevice = UseDeviceOption

!     nV = size ( S_E )
    
!     if ( UseDevice ) then
!       !$OMP OMP_TARGET_DIRECTIVE parallel do &
!       !$OMP schedule ( OMP_SCHEDULE_TARGET )
!       do iV = 1, nV
!         if ( ProperCell ( iV ) ) then      

!           S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) )  !&
! !                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
!           S_S_1 ( iV )  =  ( Xi_H ( iV )  &
!                              -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) )!&
! !                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!           S_S_2 ( iV )  =  ( Xi_H ( iV )  &
!                              -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) )!&
! !                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!           S_S_3 ( iV )  =  ( Xi_H ( iV )  &
!                              -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) )!&
! !                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )

!           if ( J_RD ( iV )  <  1.0e-3_KDR .and. SF_RD ( iV )  <  1.0e-6_KDR ) &
!           then

!             DI ( iV )  =  1.0_KDR

!             ! S_E ( iV )  =  ( J_Eq ( iV )  -  J ( iV ) )  /  dT

!             ! S_S_1 ( iV )  &
!             !   =  ( S_S_1_D ( iV )  /  Chi_H ( iV )  -  S_1 ( iV ) )  /  dT  &
!             !      -  S_S_1_D ( iV )
!             ! S_S_2 ( iV )  &
!             !   =  ( S_S_2_D ( iV )  /  Chi_H ( iV )  -  S_2 ( iV ) )  /  dT  &
!             !      -  S_S_2_D ( iV )
!             ! S_S_3 ( iV )  &
!             !   =  ( S_S_3_D ( iV )  /  Chi_H ( iV )  -  S_3 ( iV ) )  /  dT  &
!             !      -  S_S_3_D ( iV )

!           else
!             DI ( iV )  =  0.0_KDR
!           end if

!         else
!           S_E   ( iV )  =  0.0_KDR
!           S_S_1 ( iV )  =  0.0_KDR
!           S_S_2 ( iV )  =  0.0_KDR
!           S_S_3 ( iV )  =  0.0_KDR
!         end if !-- ProperCell
!       end do !-- iV
!       !$OMP end OMP_TARGET_DIRECTIVE parallel do
!     else
!       !$OMP parallel do &
!       !$OMP schedule ( OMP_SCHEDULE_HOST )
!       do iV = 1, nV
!         if ( ProperCell ( iV ) ) then      

!           S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) )  !&
! !                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
!           S_S_1 ( iV )  =  ( Xi_H ( iV )  &
!                              -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) )!&
! !                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!           S_S_2 ( iV )  =  ( Xi_H ( iV )  &
!                              -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) )!&
! !                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
!           S_S_3 ( iV )  =  ( Xi_H ( iV )  &
!                              -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) )!&
! !                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )

!           if ( J_RD ( iV )  <  1.0e-3_KDR .and. SF_RD ( iV )  <  1.0e-6_KDR ) &
!           then

!             DI ( iV )  =  1.0_KDR

!             ! S_E ( iV )  =  ( J_Eq ( iV )  -  J ( iV ) )  /  dT

!             ! S_S_1 ( iV )  &
!             !   =  ( S_S_1_D ( iV )  /  Chi_H ( iV )  -  S_1 ( iV ) )  /  dT  &
!             !      -  S_S_1_D ( iV )
!             ! S_S_2 ( iV )  &
!             !   =  ( S_S_2_D ( iV )  /  Chi_H ( iV )  -  S_2 ( iV ) )  /  dT  &
!             !      -  S_S_2_D ( iV )
!             ! S_S_3 ( iV )  &
!             !   =  ( S_S_3_D ( iV )  /  Chi_H ( iV )  -  S_3 ( iV ) )  /  dT  &
!             !      -  S_S_3_D ( iV )

!           else
!             DI ( iV )  =  0.0_KDR
!           end if

!         else
!           S_E   ( iV )  =  0.0_KDR
!           S_S_1 ( iV )  =  0.0_KDR
!           S_S_2 ( iV )  =  0.0_KDR
!           S_S_3 ( iV )  =  0.0_KDR
!         end if !-- ProperCell
!       end do !-- iV
!       !$OMP end parallel do
!     end if

!   end procedure ComputeKernel


  module procedure ComputeSimpleKernel

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
          S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) )  !&
!                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
          S_S_1 ( iV )  =  ( Xi_H ( iV )  &
                             -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) )!&
!                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
          S_S_2 ( iV )  =  ( Xi_H ( iV )  &
                             -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) )!&
!                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
          S_S_3 ( iV )  =  ( Xi_H ( iV )  &
                             -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) )!&
!                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
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
          S_E   ( iV )  =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) )!  &
!                           /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )
          S_S_1 ( iV )  =  ( Xi_H ( iV )  &
                             -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) )!&
!                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
          S_S_2 ( iV )  =  ( Xi_H ( iV )  &
                             -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) )!&
!                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
          S_S_3 ( iV )  =  ( Xi_H ( iV )  &
                             -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) )!&
!                           /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
        else
          S_E   ( iV )  =  0.0_KDR
          S_S_1 ( iV )  =  0.0_KDR
          S_S_2 ( iV )  =  0.0_KDR
          S_S_3 ( iV )  =  0.0_KDR
        end if !-- ProperCell
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeSimpleKernel


end submodule Slope_RM_I_I__Kernel
