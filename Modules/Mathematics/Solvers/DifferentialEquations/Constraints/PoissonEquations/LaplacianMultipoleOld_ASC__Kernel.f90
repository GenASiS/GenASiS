#include "Preprocessor"

submodule ( LaplacianMultipoleOld_ASC__Form ) LaplacianMultipoleOld_ASC__Kernel

  use Basics
  use LaplacianMultipoleOld_Template

  implicit none

contains


  module procedure ComputeMomentsLocalOld_CSL_Kernel

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iE, &  !-- iEquation
      iR     !-- iRadius
    real ( KDR ) :: &
      R  !-- Radius
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
!--      !$OMP  num_teams ( 4 ) thread_limit ( 16 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iC, iE, iR, R, SH_RC, SH_IC, SH_RS, SH_IS ) &
      !$OMP& reduction ( + : MyM_RC ) &
      !$OMP& reduction ( + : MyM_IC ) &
      !$OMP& reduction ( + : MyM_RS ) &
      !$OMP& reduction ( + : MyM_IS ) 
      do iC = 1, nCells

        if ( .not. IsProperCell ( iC ) ) &
          cycle

        call ComputeSolidHarmonicsKernel &
               ( CoordinateSystem, &
                 [ Center_1 ( iC ), Center_2 ( iC ), Center_3 ( iC ) ], &
                 Origin, RadialEdge, MaxDegree, nDimensions, GridError, &
                 SH_RC, SH_IC, SH_RS, SH_IS, R, iR )

        call ComputeMomentContributionsOldKernel &
               ( MyM_RC, MyM_IC, MyM_RS, MyM_IS, SH_RC, SH_IC, SH_RS, SH_IS, &
                 Source ( iC, : ), Volume ( iC ), iaSource, MaxOrder, &
                 nEquations, nAngularMomentCells, iR )

      end do !-- iC
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iC, iR, R, SH_RC, SH_IC, SH_RS, SH_IS ) &
      !$OMP& reduction ( + : MyM_RC ) &
      !$OMP& reduction ( + : MyM_IC ) &
      !$OMP& reduction ( + : MyM_RS ) &
      !$OMP& reduction ( + : MyM_IS ) 
      do iC = 1, nCells

        if ( .not. IsProperCell ( iC ) ) &
          cycle

        call ComputeSolidHarmonicsKernel &
               ( CoordinateSystem, &
                 [ Center_1 ( iC ), Center_2 ( iC ), Center_3 ( iC ) ], &
                 Origin, RadialEdge, MaxDegree, nDimensions, GridError, &
                 SH_RC, SH_IC, SH_RS, SH_IS, R, iR )

        call ComputeMomentContributionsOldKernel &
               ( MyM_RC, MyM_IC, MyM_RS, MyM_IS, SH_RC, SH_IC, SH_RS, SH_IS, &
                 Source ( iC, : ), Volume ( iC ), iaSource, MaxOrder, &
                 nEquations, nAngularMomentCells, iR )

      end do !-- iC
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure ComputeMomentsLocalOld_CSL_Kernel


end submodule LaplacianMultipoleOld_ASC__Kernel
