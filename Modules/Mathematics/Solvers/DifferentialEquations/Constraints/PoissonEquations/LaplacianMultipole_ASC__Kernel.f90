#include "Preprocessor"

submodule ( LaplacianMultipole_ASC__Form ) LaplacianMultipole_ASC__Kernel

  use Basics
  use LaplacianMultipole_Template

  implicit none

contains


  module procedure ComputeMomentsLocal_CSL_Kernel

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iR     !-- iRadius
    real ( KDR ) :: &
      R  !-- Radius

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
             ( CoordinateSystem, Center ( iC, : ), Origin, RadialEdge, &
               MaxDegree, nDimensions, GridError, SH_RC, SH_IC, SH_RS, SH_IS, &
               R, iR )

      call ComputeMomentContributionsKernel &
             ( MyM_RC, MyM_IC, MyM_RS, MyM_IS, SH_RC, SH_IC, SH_RS, SH_IS, &
               Source ( iC, : ), Volume ( iC ), iaSource, MaxOrder, &
               nEquations, nAngularMomentCells, iR )

    end do !-- iC
    !$OMP  end parallel do

  end procedure ComputeMomentsLocal_CSL_Kernel


end submodule LaplacianMultipole_ASC__Kernel
