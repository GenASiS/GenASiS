#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics
  use LaplacianMultipole_Template

  implicit none

contains


  module procedure ComputeRectangularCoordinates_CSL_Kernel

    integer ( KDI ) :: &
      iC
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    CoordinateError = .false.

    if ( UseDevice ) then
     
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC )
        do iC = 1, nC
          if ( .not. IsProperCell ( iC ) ) &
            cycle
          X ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  cos ( X_3 ( iC ) )
          Y ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  sin ( X_3 ( iC ) )
          Z ( iC )  =  X_1 ( iC )  *  cos ( X_2 ( iC ) )
        end do !-- iC
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- Host

      select case ( COORDINATE_SYSTEM )
      case ( SPHERICAL )

        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC )
        do iC = 1, nC
          if ( .not. IsProperCell ( iC ) ) &
            cycle
          X ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  cos ( X_3 ( iC ) )
          Y ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  sin ( X_3 ( iC ) )
          Z ( iC )  =  X_1 ( iC )  *  cos ( X_2 ( iC ) )
        end do !-- iC
        !$OMP  end parallel do
        
      case default
        CoordinateError = .true.
        return
      end select !-- COORDINATE_SYSTEM
      
    end if !-- UseDevice

  end procedure ComputeRectangularCoordinates_CSL_Kernel


  module procedure ComputeMomentContributions_CSL_Kernel

    integer ( KDI ) :: &
      iA, &  !-- iAngularMoment
      iM, &  !-- iOrder
      iL     !-- iDegree
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
     
      !$OMP target
 
      !$OMP end target

    else

      iA  =  0
      do iM  =  0, M
        do iL  =  iM, L

          iA  =  iA + 1

          MyM_RC ( iA, :, : )  =  1.0_KDR * iA
          MyM_IC ( iA, :, : )  =  2.0_KDR * iA
          MyM_RS ( iA, :, : )  =  3.0_KDR * iA
          MyM_IS ( iA, :, : )  =  4.0_KDR * iA

        end do !-- iL
      end do !-- iM

    end if

  end procedure ComputeMomentContributions_CSL_Kernel


end submodule LaplacianMultipole_ASC__Kernel
