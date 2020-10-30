#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics

  implicit none

contains


  module procedure ComputeMomentsLocal_CSL_S_Kernel

    !-- ComputeMomentsLocal_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iAM, &  !-- iAngularMoment
      iE     !-- iEquation
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

    else  !-- use host

      !$OMP  parallel do collapse ( 5 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iT, iP, iR, iAM, iE ) &
      !$OMP& reduction ( + : MySM )
      do iE  =  1, nE
        do iAM  =  1, nAM
          do iR  =  1,  nC ( 1 )
            do iP  =  1,  nC ( 3 )
              do iT  =  1,  nC ( 2 )

                MySM ( oR + iR, iAM, iE )  &
                  =  MySM ( oR + iR, iAM, iE )  &
                     +  dSA ( iT, iP )  * AF ( iT, iP, iAM )  &
                        *  S ( oC ( 1 )  +  iR, &
                               oC ( 2 )  +  iT, &
                               oC ( 3 )  +  iP, &
                               iE )

              end do !-- iT
            end do !-- iP
          end do !-- iR
        end do !-- iAM
      end do !-- iE
      !$OMP  end parallel do      

    end if  !-- UseDevice

  end procedure ComputeMomentsLocal_CSL_S_Kernel


end submodule LaplacianMultipole_ASC__Kernel
