#include "Preprocessor"

submodule ( Poisson_ASCG__Form ) Poisson_ASCG__Kernel

  use Basics
  
  implicit none

contains


  module procedure CombineMoments_CGS_S_Kernel

    !-- CombineMoments_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iAM, &  !-- iAngularMoment
      iE     !-- iEquation
    real ( KDR ) :: &
      RM_R_C, RM_I_C, &  !-- RadialMoment_[Regular,Irregular]_Center
      SE
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
      
      !-- The iAM loop at least must be separated to avoid an OMP reduction
      do iE  =  1, nE
        do iAM  =  1, nAM

          !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
          !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
          !$OMP private ( RM_R_C, RM_I_C ) &
          !$OMP firstprivate ( iE, iAM )
          do iP  =  1,  nC ( 3 )
            do iT  =  1,  nC ( 2 )
              do iR  =  1,  nC ( 1 )

                RM_R_C  =  0.5_KDR  &
                             *  (    RM_R ( oR + iR,     iAM, iE )  &
                                  +  RM_R ( oR + iR + 1, iAM, iE )  )
                RM_I_C  =  0.5_KDR  &
                             *  (    RM_I ( oR + iR,     iAM, iE )  &
                                  +  RM_I ( oR + iR + 1, iAM, iE )  )

                S ( oC ( 1 ) + iR, oC ( 2 ) + iT, oC ( 3 ) + iP, iE )  &
                  =  S ( oC ( 1 ) + iR, oC ( 2 ) + iT, oC ( 3 ) + iP, iE )  &
                     -  DF ( iAM )  *  AF ( iT, iP, iAM )  &
                        *  (    RF_I ( oR + iR, iAM )  *  RM_R_C &
                             +  RF_R ( oR + iR, iAM )  *  RM_I_C )

              end do !-- iR
            end do !-- iT
          end do !-- iP
          !$OMP end OMP_TARGET_DIRECTIVE parallel do      

        end do !-- iAM
      end do !-- iE
      
    else  !-- use host

      !-- The iAM loop at least must be separated to avoid an OMP reduction
      do iE  =  1, nE
        do iAM  =  1, nAM

          !$OMP parallel do collapse ( 3 ) &
          !$OMP schedule ( OMP_SCHEDULE_HOST ) &
          !$OMP private ( RM_R_C, RM_I_C ) &
          !$OMP firstprivate ( iE, iAM )
          do iP  =  1,  nC ( 3 )
            do iT  =  1,  nC ( 2 )
              do iR  =  1,  nC ( 1 )

                RM_R_C  =  0.5_KDR  &
                             *  (    RM_R ( oR + iR,     iAM, iE )  &
                                  +  RM_R ( oR + iR + 1, iAM, iE )  )
                RM_I_C  =  0.5_KDR  &
                             *  (    RM_I ( oR + iR,     iAM, iE )  &
                                  +  RM_I ( oR + iR + 1, iAM, iE )  )

                S ( oC ( 1 ) + iR, oC ( 2 ) + iT, oC ( 3 ) + iP, iE )  &
                  =  S ( oC ( 1 ) + iR, oC ( 2 ) + iT, oC ( 3 ) + iP, iE )  &
                     -  DF ( iAM )  *  AF ( iT, iP, iAM )  &
                        *  (    RF_I ( oR + iR, iAM )  *  RM_R_C &
                             +  RF_R ( oR + iR, iAM )  *  RM_I_C )

              end do !-- iR
            end do !-- iT
          end do !-- iP
          !$OMP end parallel do      

        end do !-- iAM
      end do !-- iE

    end if  !-- UseDevice

  end procedure CombineMoments_CGS_S_Kernel


end submodule Poisson_ASCG__Kernel
