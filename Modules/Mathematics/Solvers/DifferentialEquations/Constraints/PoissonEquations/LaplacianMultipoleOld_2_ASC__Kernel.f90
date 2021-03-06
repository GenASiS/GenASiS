#include "Preprocessor"

submodule ( LaplacianMultipoleOld_2_ASC__Form ) LaplacianMultipoleOld_2_ASC__Kernel

  use Basics

  implicit none

contains


  module procedure Compute_RC_CSL_S_Kernel

    !-- Compute_RectangularCoordinates_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iC
    logical ( KDL ) :: &
      UseDevice

    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
     
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC )
      do iC  =  1, nC

        if ( IsProperCell ( iC ) ) then
          X ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  &
                                   *  cos ( X_3 ( iC ) )
          Y ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  &
                                   *  sin ( X_3 ( iC ) )
          Z ( iC )  =  X_1 ( iC )  *  cos ( X_2 ( iC ) )
        else
          X ( iC )  =  0.0_KDR
          Y ( iC )  =  0.0_KDR
          Z ( iC )  =  0.0_KDR
        end if

        D   ( iC )  =  X_1 ( iC )
        D_2 ( iC )  =  X_1 ( iC ) ** 2

      end do !-- iC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC )
      do iC  =  1, nC

        if ( IsProperCell ( iC ) ) then
          X ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  &
                                   *  cos ( X_3 ( iC ) )
          Y ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  &
                                   *  sin ( X_3 ( iC ) )
          Z ( iC )  =  X_1 ( iC )  *  cos ( X_2 ( iC ) )
        else
          X ( iC )  =  0.0_KDR
          Y ( iC )  =  0.0_KDR
          Z ( iC )  =  0.0_KDR
        end if

        D   ( iC )  =  X_1 ( iC )
        D_2 ( iC )  =  X_1 ( iC ) ** 2

      end do !-- iC
      !$OMP  end parallel do
        
    end if !-- UseDevice

  end procedure Compute_RC_CSL_S_Kernel


  module procedure Compute_SH_0_0_CSL_Kernel

    !-- Compute_SolidHarmonics_0_0_ChartSingleLevel_Kernel
    !-- See Flash 4 User's Guide p. 129

    !-- ( L, M ) = ( iM, iM )
    !-- Note iL = iM

    integer ( KDI ) :: &
      iC, jC, kC  !-- iCell, etc.
    logical ( KDL ) :: &
      UseDevice

    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC, jC, kC ) 
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR
            I_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR / sqrt ( D_2 ( iC, jC, kC ) )
            R_S ( iC, jC, kC, iSH_0 )  =  0.0_KDR
            I_S ( iC, jC, kC, iSH_0 )  =  0.0_KDR

            R_C ( iC, jC, kC, iSH_PD )  =  R_C ( iC, jC, kC, iSH_0 )
            I_C ( iC, jC, kC, iSH_PD )  =  I_C ( iC, jC, kC, iSH_0 )
            R_S ( iC, jC, kC, iSH_PD )  =  R_S ( iC, jC, kC, iSH_0 )
            I_S ( iC, jC, kC, iSH_PD )  =  I_S ( iC, jC, kC, iSH_0 )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR
            I_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR / sqrt ( D_2 ( iC, jC, kC ) )
            R_S ( iC, jC, kC, iSH_0 )  =  0.0_KDR
            I_S ( iC, jC, kC, iSH_0 )  =  0.0_KDR

            R_C ( iC, jC, kC, iSH_PD )  =  R_C ( iC, jC, kC, iSH_0 )
            I_C ( iC, jC, kC, iSH_PD )  =  I_C ( iC, jC, kC, iSH_0 )
            R_S ( iC, jC, kC, iSH_PD )  =  R_S ( iC, jC, kC, iSH_0 )
            I_S ( iC, jC, kC, iSH_PD )  =  I_S ( iC, jC, kC, iSH_0 )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure Compute_SH_0_0_CSL_Kernel


  module procedure Compute_SH_iM_iM_CSL_Kernel

    !-- Compute_SolidHarmonics_M_M_ChartSingleLevel_Kernel
    !-- See Flash 4 User's Guide p. 129

    !-- ( L, M ) = ( iM, iM )
    !-- Note iL = iM

    integer ( KDI ) :: &
      iC, jC, kC  !-- iCell, etc.
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  &
              =  - (    X ( iC, jC, kC )  *  R_C ( iC, jC, kC, iSH_PD )  &
                     -  Y ( iC, jC, kC )  *  R_S ( iC, jC, kC, iSH_PD ) )  &
                   /  ( 2 * iM )

            I_C ( iC, jC, kC, iSH_0 )  &
              =  - ( 2 * iM - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  (    X ( iC, jC, kC )  *  I_C ( iC, jC, kC, iSH_PD )  &
                        -  Y ( iC, jC, kC )  *  I_S ( iC, jC, kC, iSH_PD ) )

            R_S ( iC, jC, kC, iSH_0 )  &
              =  - (    Y ( iC, jC, kC )  *  R_C ( iC, jC, kC, iSH_PD )  &
                     +  X ( iC, jC, kC )  *  R_S ( iC, jC, kC, iSH_PD ) )  &
                   /  ( 2 * iM )

            I_S ( iC, jC, kC, iSH_0 )  &
              =  - ( 2 * iM - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  (    Y ( iC, jC, kC )  *  I_C ( iC, jC, kC, iSH_PD )  &
                        +  X ( iC, jC, kC )  *  I_S ( iC, jC, kC, iSH_PD ) )

            R_C ( iC, jC, kC, iSH_PD )  =  R_C ( iC, jC, kC, iSH_0 )
            I_C ( iC, jC, kC, iSH_PD )  =  I_C ( iC, jC, kC, iSH_0 )
            R_S ( iC, jC, kC, iSH_PD )  =  R_S ( iC, jC, kC, iSH_0 )
            I_S ( iC, jC, kC, iSH_PD )  =  I_S ( iC, jC, kC, iSH_0 )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  &
              =  - (    X ( iC, jC, kC )  *  R_C ( iC, jC, kC, iSH_PD )  &
                     -  Y ( iC, jC, kC )  *  R_S ( iC, jC, kC, iSH_PD ) )  &
                   /  ( 2 * iM )

            I_C ( iC, jC, kC, iSH_0 )  &
              =  - ( 2 * iM - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  (    X ( iC, jC, kC )  *  I_C ( iC, jC, kC, iSH_PD )  &
                        -  Y ( iC, jC, kC )  *  I_S ( iC, jC, kC, iSH_PD ) )

            R_S ( iC, jC, kC, iSH_0 )  &
              =  - (    Y ( iC, jC, kC )  *  R_C ( iC, jC, kC, iSH_PD )  &
                     +  X ( iC, jC, kC )  *  R_S ( iC, jC, kC, iSH_PD ) )  &
                   /  ( 2 * iM )

            I_S ( iC, jC, kC, iSH_0 )  &
              =  - ( 2 * iM - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  (    Y ( iC, jC, kC )  *  I_C ( iC, jC, kC, iSH_PD )  &
                        +  X ( iC, jC, kC )  *  I_S ( iC, jC, kC, iSH_PD ) )

            R_C ( iC, jC, kC, iSH_PD )  =  R_C ( iC, jC, kC, iSH_0 )
            I_C ( iC, jC, kC, iSH_PD )  =  I_C ( iC, jC, kC, iSH_0 )
            R_S ( iC, jC, kC, iSH_PD )  =  R_S ( iC, jC, kC, iSH_0 )
            I_S ( iC, jC, kC, iSH_PD )  =  I_S ( iC, jC, kC, iSH_0 )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure Compute_SH_iM_iM_CSL_Kernel


  module procedure Compute_SH_iL_iM_1_CSL_Kernel

    !-- Compute_SolidHarmonics_L_M_1_ChartSingleLevel_Kernel
    !-- See Flash 4 User's Guide p. 129

    !-- ( L, M ) = ( iM + 1, iM )
    !-- Note iL = iM + 1

    integer ( KDI ) :: &
      iC, jC, kC  !-- iCell, etc.
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  &
              =  Z ( iC, jC, kC )  *  R_C ( iC, jC, kC, iSH_1 )

            I_C ( iC, jC, kC, iSH_0 )  &
              =  ( 2 * ( iM + 1 ) - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  Z ( iC, jC, kC )  *  I_C ( iC, jC, kC, iSH_1 )

            R_S ( iC, jC, kC, iSH_0 )  &
              =  Z ( iC, jC, kC )  *  R_S ( iC, jC, kC, iSH_1 )

            I_S ( iC, jC, kC, iSH_0 )  &
              =  ( 2 * ( iM + 1 ) - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  Z ( iC, jC, kC )  *  I_S ( iC, jC, kC, iSH_1 )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  &
              =  Z ( iC, jC, kC )  *  R_C ( iC, jC, kC, iSH_1 )

            I_C ( iC, jC, kC, iSH_0 )  &
              =  ( 2 * ( iM + 1 ) - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  Z ( iC, jC, kC )  *  I_C ( iC, jC, kC, iSH_1 )

            R_S ( iC, jC, kC, iSH_0 )  &
              =  Z ( iC, jC, kC )  *  R_S ( iC, jC, kC, iSH_1 )

            I_S ( iC, jC, kC, iSH_0 )  &
              =  ( 2 * ( iM + 1 ) - 1 )  /  D_2 ( iC, jC, kC )  &
                   *  Z ( iC, jC, kC )  *  I_S ( iC, jC, kC, iSH_1 )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure Compute_SH_iL_iM_1_CSL_Kernel


  module procedure Compute_SH_iL_iM_2_CSL_Kernel

    !-- Compute_SolidHarmonics_L_M_2_ChartSingleLevel_Kernel
    !-- See Flash 4 User's Guide p. 129

    !-- ( L, M ) = ( iL, iM )

    integer ( KDI ) :: &
      iC, jC, kC  !-- iCell, etc.
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  R_C ( iC, jC, kC, iSH_1 )  &
                   -  D_2 ( iC, jC, kC )  &
                      *  R_C ( iC, jC, kC, iSH_2 ) )  &
                 /  ( ( iL + iM ) * ( iL - iM ) )

            I_C ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  I_C ( iC, jC, kC, iSH_1 )  &
                   -  ( ( iL - 1 ) ** 2 - iM ** 2 )  &
                      *  I_C ( iC, jC, kC, iSH_2 ) )  &
                 /  D_2 ( iC, jC, kC )

            R_S ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  R_S ( iC, jC, kC, iSH_1 )  &
                   -  D_2 ( iC, jC, kC )  &
                      *  R_S ( iC, jC, kC, iSH_2 ) )  &
                 /  ( ( iL + iM ) * ( iL - iM ) )

            I_S ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  I_S ( iC, jC, kC, iSH_1 )  &
                   -  ( ( iL - 1 ) ** 2 - iM ** 2 )  &
                      *  I_S ( iC, jC, kC, iSH_2 ) )  &
                 /  D_2 ( iC, jC, kC )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
        do jC  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
          do iC  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

            R_C ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  R_C ( iC, jC, kC, iSH_1 )  &
                   -  D_2 ( iC, jC, kC )  &
                      *  R_C ( iC, jC, kC, iSH_2 ) )  &
                 /  ( ( iL + iM ) * ( iL - iM ) )

            I_C ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  I_C ( iC, jC, kC, iSH_1 )  &
                   -  ( ( iL - 1 ) ** 2 - iM ** 2 )  &
                      *  I_C ( iC, jC, kC, iSH_2 ) )  &
                 /  D_2 ( iC, jC, kC )

            R_S ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  R_S ( iC, jC, kC, iSH_1 )  &
                   -  D_2 ( iC, jC, kC )  &
                      *  R_S ( iC, jC, kC, iSH_2 ) )  &
                 /  ( ( iL + iM ) * ( iL - iM ) )

            I_S ( iC, jC, kC, iSH_0 )  &
              =  ( ( 2 * iL - 1 )  *  Z ( iC, jC, kC )  &
                      *  I_S ( iC, jC, kC, iSH_1 )  &
                   -  ( ( iL - 1 ) ** 2 - iM ** 2 )  &
                      *  I_S ( iC, jC, kC, iSH_2 ) )  &
                 /  D_2 ( iC, jC, kC )

          end do !-- iC
        end do !-- jC
      end do !-- kC
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure Compute_SH_iL_iM_2_CSL_Kernel


  module procedure ComputeMomentLocal_CSL_S_Kernel

    !-- ComputeMomentLocal_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iE  !-- iEquation
    real ( KDR ) :: &
      M_RC, M_IC, M_RS, M_IS
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DISTRIBUTE_DIRECTIVE collapse ( 2 ) &
      !$OMP& OMP_TARGET_DISTRIBUTE_SCHEDULE &  
      !$OMP& private ( iR, iE, M_RC, M_IC, M_RS, M_IS )
      do iE  =  1, nE
        do iR  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

          M_RC  =  0.0_KDR
          M_IC  =  0.0_KDR
          M_RS  =  0.0_KDR
          M_IS  =  0.0_KDR

          !$OMP  parallel do collapse ( 2 ) &
          !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iT, iP ) &
          !$OMP& reduction ( + : M_RC ) &
          !$OMP& reduction ( + : M_IC ) &
          !$OMP& reduction ( + : M_RS ) &
          !$OMP& reduction ( + : M_IS )
          do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
            do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )

              M_RC  =  M_RC  &
                       +  SH_RC ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                          *  dV ( iR, iT, iP )

              M_IC  =  M_IC  &
                       +  SH_IC ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                          *  dV ( iR, iT, iP )

              M_RS  =  M_RS  &
                       +  SH_RS ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                          *  dV ( iR, iT, iP )

              M_IS  =  M_IS  &
                       +  SH_IS ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                          *  dV ( iR, iT, iP )

            end do !-- iT
          end do !-- iP
          !$OMP  end parallel do

          MyM_RC ( oR - oC ( 1 ) + iR, iE )  =  M_RC
          MyM_IC ( oR - oC ( 1 ) + iR, iE )  =  M_IC
          MyM_RS ( oR - oC ( 1 ) + iR, iE )  =  M_RS
          MyM_IS ( oR - oC ( 1 ) + iR, iE )  =  M_IS

        end do !-- iR
      end do !-- iE
      !$OMP  end OMP_TARGET_DISTRIBUTE_DIRECTIVE  

    else !-- use host

      !$OMP  parallel do collapse ( 4 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iR, iT, iP, iE ) &
      !$OMP& reduction ( + : MyM_RC ) &
      !$OMP& reduction ( + : MyM_IC ) &
      !$OMP& reduction ( + : MyM_RS ) &
      !$OMP& reduction ( + : MyM_IS )
      do iE  =  1, nE
        do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
          do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
            do iR  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

              MyM_RC ( oR - oC ( 1 ) + iR, iE )  &
                =  MyM_RC ( oR - oC ( 1 ) + iR, iE )  &
                   +  SH_RC ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_IC ( oR - oC ( 1 ) + iR, iE )  &
                =  MyM_IC ( oR - oC ( 1 ) + iR, iE )  &
                   +  SH_IC ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_RS ( oR - oC ( 1 ) + iR, iE )  &
                =  MyM_RS ( oR - oC ( 1 ) + iR, iE )  &
                   +  SH_RS ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_IS ( oR - oC ( 1 ) + iR, iE )  &
                =  MyM_IS ( oR - oC ( 1 ) + iR, iE )  &
                   +  SH_IS ( iR, iT, iP )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

            end do !-- iR
          end do !-- iT
        end do !-- iP
      end do !-- iE
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure ComputeMomentLocal_CSL_S_Kernel


end submodule LaplacianMultipoleOld_2_ASC__Kernel
