#include "Preprocessor"

submodule ( Laplacian_M_H__Form ) Laplacian_M_H__Kernel

  use Basics

  implicit none

contains


  module procedure ComputeRadialMomentsKernel

    integer ( KDI ) :: &
      iR, &  !-- iRadius
      iAM, &  !-- iAngularMoment
      iE     !-- iEquation
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    !-- No OMP directive on iR loops because of dependence.

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iR )
      do iE  =  1, nE
        do iAM  =  1, nAM

          RM_R ( 1, iAM, iE )  =  0.0_KDR
          do iR  =  1, nR
            RM_R ( iR + 1, iAM, iE )  &
              =  RM_R ( iR, iAM, iE )  &
                 +  dR33 ( iR )  *  RF_R ( iR, iAM )  *  AM ( iR, iAM, iE )
          end do !-- iR

          RM_I ( nR + 1, iAM, iE )  =  0.0_KDR
          do iR  =  nR, 1, -1
            RM_I ( iR, iAM, iE )  &
              =  RM_I ( iR + 1, iAM, iE )  &
                 +  dR33 ( iR )  *  RF_I ( iR, iAM )  *  AM ( iR, iAM, iE )
          end do !-- iR

        end do !-- iAM
      end do !-- iE
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else  !-- use host

      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iR )
      do iE  =  1, nE
        do iAM  =  1, nAM

          RM_R ( 1, iAM, iE )  =  0.0_KDR
          do iR  =  1, nR
            RM_R ( iR + 1, iAM, iE )  &
              =  RM_R ( iR, iAM, iE )  &
                 +  dR33 ( iR )  *  RF_R ( iR, iAM )  *  AM ( iR, iAM, iE )
          end do !-- iR

          RM_I ( nR + 1, iAM, iE )  =  0.0_KDR
          do iR  =  nR, 1, -1
            RM_I ( iR, iAM, iE )  &
              =  RM_I ( iR + 1, iAM, iE )  &
                 +  dR33 ( iR )  *  RF_I ( iR, iAM )  *  AM ( iR, iAM, iE )
          end do !-- iR

        end do !-- iAM
      end do !-- iE
      !$OMP  end parallel do      

    end if  !-- UseDevice

  end procedure ComputeRadialMomentsKernel


end submodule Laplacian_M_H__Kernel
