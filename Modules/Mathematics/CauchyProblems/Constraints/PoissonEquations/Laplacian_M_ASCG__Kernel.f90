#include "Preprocessor"

submodule ( Laplacian_M_ASCG__Form ) Laplacian_M_ASCG__Kernel

  use Basics

  implicit none

contains


  module procedure ComputeAngularMomentsLocal_CGS_S_Kernel

    !-- ComputeAngularMomentsLocal_ChartGridStructured_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iAM, &  !-- iAngularMoment
      iE     !-- iEquation
    real ( KDR ) :: &
      MyAME  !-- MyAngularMomentElement
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DISTRIBUTE_DIRECTIVE collapse ( 3 ) &
      !$OMP OMP_TARGET_DISTRIBUTE_SCHEDULE private ( MyAME )
      do iE  =  1, nE
        do iAM  =  1, nAM
          do iR  =  1,  nC ( 1 )
            
            MyAME = 0.0_KDR
             
            !$OMP parallel do collapse ( 2 ) &
            !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iT, iP ) &
            !$OMP firstprivate ( iR, iAM, iE ) &
            !$OMP reduction ( + : MyAME )
            do iP  =  1,  nC ( 3 )
              do iT  =  1,  nC ( 2 )

                MyAME  =  MyAME &
                          +  dSA ( iT, iP )  * AF ( iT, iP, iAM )  &
                             *  S ( oC ( 1 )  +  iR, &
                                    oC ( 2 )  +  iT, &
                                    oC ( 3 )  +  iP, &
                                    iE )

              end do !-- iT
            end do !-- iP
            !$OMP end parallel do

            MyAM ( oR + iR, iAM, iE )  =  MyAME

          end do !-- iR
        end do !-- iAM
      end do !-- iE
      !$OMP end OMP_TARGET_DISTRIBUTE_DIRECTIVE

    else  !-- use host

      !$OMP parallel do collapse ( 5 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP reduction ( + : MyAM )
      do iE  =  1, nE
        do iAM  =  1, nAM
          do iP  =  1,  nC ( 3 )
            do iT  =  1,  nC ( 2 )
              do iR  =  1,  nC ( 1 )

                MyAM ( oR + iR, iAM, iE )  &
                  =  MyAM ( oR + iR, iAM, iE )  &
                     +  dSA ( iT, iP )  * AF ( iT, iP, iAM )  &
                        *  S ( oC ( 1 )  +  iR, &
                               oC ( 2 )  +  iT, &
                               oC ( 3 )  +  iP, &
                               iE )

              end do !-- iR
            end do !-- iT
          end do !-- iP
        end do !-- iAM
      end do !-- iE
      !$OMP  end parallel do      

    end if  !-- UseDevice

  end procedure ComputeAngularMomentsLocal_CGS_S_Kernel


end submodule Laplacian_M_ASCG__Kernel
