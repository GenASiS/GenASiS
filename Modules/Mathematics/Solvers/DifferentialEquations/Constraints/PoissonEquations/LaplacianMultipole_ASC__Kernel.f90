#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics
  use LaplacianMultipole_Template

  implicit none

contains


  module procedure ComputeRectangularCoordinates_CSL_Kernel

    integer ( KDI ) :: &
      iC
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    CoordinateError = .false.

    if ( UseDevice ) then
     
      select case ( COORDINATE_SYSTEM )
      case ( SPHERICAL )

        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC )
        do iC  =  1, nC
          if ( .not. IsProperCell ( iC ) ) &
            cycle
          X ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  cos ( X_3 ( iC ) )
          Y ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  sin ( X_3 ( iC ) )
          Z ( iC )  =  X_1 ( iC )  *  cos ( X_2 ( iC ) )
        end do !-- iC
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      case default
        CoordinateError = .true.
        return
      end select !-- COORDINATE_SYSTEM
      
    else !-- use host

      select case ( COORDINATE_SYSTEM )
      case ( SPHERICAL )

        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC )
        do iC  =  1, nC
          if ( .not. IsProperCell ( iC ) ) &
            cycle
          X ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  cos ( X_3 ( iC ) )
          Y ( iC )  =  X_1 ( iC )  *  sin ( X_2 ( iC ) )  *  sin ( X_3 ( iC ) )
          Z ( iC )  =  X_1 ( iC )  *  cos ( X_2 ( iC ) )
        end do !-- iC
        !$OMP  end parallel do
        
      case default
        CoordinateError = .true.
        return
      end select !-- COORDINATE_SYSTEM
      
    end if !-- UseDevice

  end procedure ComputeRectangularCoordinates_CSL_Kernel


  module procedure ComputeMomentContributions_CSL_Kernel

    integer ( KDI ) :: &
      iA, &  !-- iAngularMoment
      iM, &  !-- iOrder
      iL, &  !-- iDegree
      iSH_0, &  !-- iSolidHarmonic_Current
      iSH_1, &  !-- iSolidHarmonic_Previous_1
      iSH_2, &  !-- iSolidHarmonic_Previous_2
      iSH_PD    !-- iSolidHarmonic_PreviousDiagonal
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
     
      !$OMP target
 
      !$OMP end target

    else

          iA  =  1
       iSH_0  =  1
       iSH_1  =  2
       iSH_2  =  3
      iSH_PD  =  4

      do iM  =  0, M

        if ( iM  ==  0 ) then
          call ComputeSolidHarmonics_0_0_CSL_Kernel &
                 ( SH_RC, SH_IC, SH_RS, SH_IS, X, Y, Z, &
                   nCells, iSH_0, iSH_PD, UseDeviceOption )
        else
        end if

        do iL  =  iM, L

          ! MyM_RC ( iA, :, : )  =  1.0_KDR * iA
          ! MyM_IC ( iA, :, : )  =  2.0_KDR * iA
          ! MyM_RS ( iA, :, : )  =  3.0_KDR * iA
          ! MyM_IS ( iA, :, : )  =  4.0_KDR * iA

             iA  =  iA + 1
          iSH_0  =  mod ( iSH_0, 3 )  +  1  
          iSH_1  =  mod ( iSH_1, 3 )  +  1  
          iSH_2  =  mod ( iSH_2, 3 )  +  1  

        end do !-- iL
      end do !-- iM

    end if

  end procedure ComputeMomentContributions_CSL_Kernel


  module procedure ComputeSolidHarmonics_0_0_CSL_Kernel

    !$OMP declare target

    integer ( KDI ) :: &
      iC, jC, kC  !-- iCell, etc.
    real ( KDR ) :: &
      D_2
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iC, jC, kC, D_2 )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

            D_2  =      X ( iC, jC, kC )  *  X ( iC, jC, kC )  &
                    +   Y ( iC, jC, kC )  *  Y ( iC, jC, kC )  &
                    +   Z ( iC, jC, kC )  *  Z ( iC, jC, kC )

            R_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR
            I_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR / sqrt ( D_2 )
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
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, jC, kC, D_2 )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

            D_2  =      X ( iC, jC, kC )  *  X ( iC, jC, kC )  &
                    +   Y ( iC, jC, kC )  *  Y ( iC, jC, kC )  &
                    +   Z ( iC, jC, kC )  *  Z ( iC, jC, kC )

            R_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR
            I_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR / sqrt ( D_2 )
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

  end procedure ComputeSolidHarmonics_0_0_CSL_Kernel


  module procedure ComputeSolidHarmonics_CSL_Kernel

    !$OMP declare target

    if ( nDimensions < 3 ) then
    else !-- nDimensions == 3 
    end if !-- nDimensions

  end procedure ComputeSolidHarmonics_CSL_Kernel


  module procedure Compute_SH_CSL_C_M_0_Kernel

    !$OMP declare target

    integer ( KDI ) :: &
      iC, jC, kC, &  !-- iCell, etc.
      iCV  !-- iCellValue
    real ( KDR ) :: &
      D_2
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    iCV  =  0

    if ( UseDevice ) then

    else !-- use host

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iC, jC, kC, D_2 ) firstprivate ( iCV )
      do kC  =  1, nCells ( 3 )
        do jC  =  1, nCells ( 2 )
          do iC  =  1, nCells ( 1 )

            iCV = iCV + 1
            if ( .not. IsProperCell ( iCV ) ) &
              cycle

            D_2  =      X ( iCV )  *  X ( iCV )  &
                    +   Z ( iCV )  *  Z ( iCV )

            if ( iM  ==  0 ) then

              R_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR
              I_C ( iC, jC, kC, iSH_0 )  =  1.0_KDR / sqrt ( D_2 )

            else

              R_C ( iC, jC, kC, iSH_0 )  &
                =  - ( X ( iCV )  *  R_C ( iC, jC, kC, iSH_PD ) ) &
                     /  ( 2 * iM )

              I_C ( iC, jC, kC, iSH_0 )  &
                =  - ( 2 * iM - 1 )  &
                     *  ( X ( iCV )  *  I_C ( iC, jC, kC, iSH_PD ) )  /  D_2

            end if
            
            R_C ( iC, jC, kC, iSH_PD )  =  R_C ( iC, jC, kC, iSH_0 )
            I_C ( iC, jC, kC, iSH_PD )  =  I_C ( iC, jC, kC, iSH_0 )

            if ( iM  <  L ) then
            

            end if !-- iM < L

          end do !-- iC
        end do !-- jC
      end do !-- kC

    end if !-- UseDevice

  end procedure Compute_SH_CSL_C_M_0_Kernel


end submodule LaplacianMultipole_ASC__Kernel
