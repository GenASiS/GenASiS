#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics
  use LaplacianMultipole_Template

  implicit none

contains


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
