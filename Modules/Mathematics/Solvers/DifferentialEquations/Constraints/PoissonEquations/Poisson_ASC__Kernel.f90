#include "Preprocessor"

submodule ( Poisson_ASC__Form ) Poisson_ASC__Kernel

  use OMP_LIB
  use Basics 
  
  implicit none

contains


  module procedure SolveCells_CSL_Kernel

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iT, &  !-- iThread
      iR     !-- iRadius
    real ( KDR ) :: &
      R, &  !-- Radius
      S     !-- Solution

    call Show ( '>>> Hello world from SolveCells_CSL_Kernel' )

    !$OMP  parallel do &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iC, iT, iR, R, S )
    do iC = 1, nCells

      if ( .not. IsProperCell ( iC ) ) &
        cycle

      iT = OMP_GET_THREAD_NUM ( )

      call ComputeSolidHarmonics ( )!&
!             ( CoordinateSystem, IsProperCell, Center ( iC, : ), nDimensions, &
!               R, iR )

    end do !-- iC
    !$OMP  end parallel do

  end procedure SolveCells_CSL_Kernel 


end submodule Poisson_ASC__Kernel
