#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics

  implicit none

contains


  module procedure ComputeMomentsLocal_CSL_S_Kernel

    !-- ComputeMomentsLocal_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iM, &  !-- iMoment
      iE     !-- iEquation
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

    else  !-- use host

      ! !$OMP  parallel do collapse ( 4 ) &
      ! !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iR, iT, iP, iE ) &
      ! !$OMP& reduction ( + : MyM )
      ! do iE  =  1, nE
      !   do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
      !     do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
      !       do iR  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

      !         MyM_RC ( oR - oC ( 1 ) + iR, iE )  &
      !           =  MyM_RC ( oR - oC ( 1 ) + iR, iE )  &
      !              +  SH_RC ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
      !                 *  dV ( iR, iT, iP )

      !         MyM_IC ( oR - oC ( 1 ) + iR, iE )  &
      !           =  MyM_IC ( oR - oC ( 1 ) + iR, iE )  &
      !              +  SH_IC ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
      !                 *  dV ( iR, iT, iP )

      !         MyM_RS ( oR - oC ( 1 ) + iR, iE )  &
      !           =  MyM_RS ( oR - oC ( 1 ) + iR, iE )  &
      !              +  SH_RS ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
      !                 *  dV ( iR, iT, iP )

      !         MyM_IS ( oR - oC ( 1 ) + iR, iE )  &
      !           =  MyM_IS ( oR - oC ( 1 ) + iR, iE )  &
      !              +  SH_IS ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
      !                 *  dV ( iR, iT, iP )

      !       end do !-- iR
      !     end do !-- iT
      !   end do !-- iP
      ! end do !-- iE
      ! !$OMP  end parallel do      

    end if  !-- UseDevice

  end procedure ComputeMomentsLocal_CSL_S_Kernel


end submodule LaplacianMultipole_ASC__Kernel
