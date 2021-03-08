#include "Preprocessor"

submodule ( LaplacianMultipoleOld_2__Template ) LaplacianMultipoleOld_2__Kernel

  use Basics

  implicit none

contains


  module procedure AddMomentShellsKernel

    integer ( KDI ) :: &
      iR, &  !-- iRadius
      iE, &  !-- iEquation
      iA     !-- iAngularMoment
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    !-- No OMP directive on iR loops because of dependence.

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iR, iE, iA )
      do iA  =  1, nA
        do iE  =  1, nE

          do iR  =  2, nR

            M_RC ( iR, iE, iA )  &
              =  M_RC ( iR - 1, iE, iA )  +  M_RC ( iR, iE, iA )

            M_RS ( iR, iE, iA )  &
              =  M_RS ( iR - 1, iE, iA )  +  M_RS ( iR, iE, iA )

          end do !-- iR

          do iR  =  nR - 1, 1, -1

            M_IC ( iR, iE, iA )  &
              =  M_IC ( iR + 1, iE, iA )  +  M_IC ( iR, iE, iA )
 
            M_IS ( iR, iE, iA )  &
              =  M_IS ( iR + 1, iE, iA )  +  M_IS ( iR, iE, iA )

          end do !-- iR

        end do !-- iE
      end do !-- iA
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do collapse ( 2 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iR, iE, iA )
      do iA  =  1, nA
        do iE  =  1, nE

          do iR  =  2, nR

            M_RC ( iR, iE, iA )  &
              =  M_RC ( iR - 1, iE, iA )  +  M_RC ( iR, iE, iA )

            M_RS ( iR, iE, iA )  &
              =  M_RS ( iR - 1, iE, iA )  +  M_RS ( iR, iE, iA )

          end do !-- iR

          do iR  =  nR - 1, 1, -1

            M_IC ( iR, iE, iA )  &
              =  M_IC ( iR + 1, iE, iA )  +  M_IC ( iR, iE, iA )
 
            M_IS ( iR, iE, iA )  &
              =  M_IS ( iR + 1, iE, iA )  +  M_IS ( iR, iE, iA )

          end do !-- iR

        end do !-- iE
      end do !-- iA
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure AddMomentShellsKernel


end submodule LaplacianMultipoleOld_2__Kernel
