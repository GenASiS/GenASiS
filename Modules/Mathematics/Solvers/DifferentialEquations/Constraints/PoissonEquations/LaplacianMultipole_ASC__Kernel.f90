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
      !$OMP& reduction ( + : MyM )
      do iE  =  1, nE
        do iAM  =  1, nAM
          do iR  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )
            do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
              do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )

                MyM ( oR - oC ( 1 ) + iR, iAM, iE )  &
                  =  MyM ( oR - oC ( 1 ) + iR, iAM, iE )  &
                     +  dSA ( iT, iP )  * AF ( iT, iP, iAM )  &
                        *  S ( iR, iT, iP, iE )

              end do !-- iT
            end do !-- iP
          end do !-- iR
        end do !-- iAM
      end do !-- iE
      !$OMP  end parallel do      

    end if  !-- UseDevice

  end procedure ComputeMomentsLocal_CSL_S_Kernel


end submodule LaplacianMultipole_ASC__Kernel
