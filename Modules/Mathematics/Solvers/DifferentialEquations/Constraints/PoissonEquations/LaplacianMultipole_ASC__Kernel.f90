#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics
  use LaplacianMultipole_Template

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

        D_2 ( iC )  =  X ( iC ) ** 2  +  Y ( iC ) ** 2  +  Z ( iC ) ** 2

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

        D_2 ( iC )  =  X ( iC ) ** 2  +  Y ( iC ) ** 2  +  Z ( iC ) ** 2

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


  module procedure Sum_MC_CSL_S_Kernel

    !-- Sum_MomentContribtions_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iE  !-- iEquation
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 4 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iR, iT, iP, iE )
      do iE  =  1, nE
        do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
          do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
            do iR  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

              MyM_RC ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_RC ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_RC ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_IC ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_IC ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_IC ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_RS ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_RS ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_RS ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_IS ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_IS ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_IS ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

            end do !-- iR
          end do !-- iT
        end do !-- iP
      end do !-- iE
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do collapse ( 4 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iR, iT, iP, iE )
      do iE  =  1, nE
        do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
          do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
            do iR  =  oC ( 1 )  +  1,  oC ( 1 )  +  nC ( 1 )

              MyM_RC ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_RC ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_RC ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_IC ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_IC ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_IC ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_RS ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_RS ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_RS ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

              MyM_IS ( iA, iR - oC ( 1 ), iE )  &
                =  MyM_IS ( iA, iR - oC ( 1 ), iE )  &
                   +  SH_IS ( iR, iT, iP, iSH_0 )  *  S ( iR, iT, iP, iE )  &
                      *  dV ( iR, iT, iP )

            end do !-- iR
          end do !-- iT
        end do !-- iP
      end do !-- iE
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure Sum_MC_CSL_S_Kernel


end submodule LaplacianMultipole_ASC__Kernel
