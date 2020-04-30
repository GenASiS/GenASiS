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
    real ( KDR ) :: &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice

    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then
     
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC ) &
      !$OMP& firstprivate ( SqrtTiny )
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

        D_2 ( iC )  =  max (     X ( iC )  *  X ( iC )  &
                             +   Y ( iC )  *  Y ( iC )  &
                             +   Z ( iC )  *  Z ( iC ),  &
                             SqrtTiny )

      end do !-- iC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC ) &
      !$OMP& firstprivate ( SqrtTiny )
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

        D_2 ( iC )  =  max (     X ( iC )  *  X ( iC )  &
                             +   Y ( iC )  *  Y ( iC )  &
                             +   Z ( iC )  *  Z ( iC ),  &
                             SqrtTiny )

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
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

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


  module procedure SumMomentContributions_CSL_SphericalKernel

  end procedure SumMomentContributions_CSL_SphericalKernel


end submodule LaplacianMultipole_ASC__Kernel
