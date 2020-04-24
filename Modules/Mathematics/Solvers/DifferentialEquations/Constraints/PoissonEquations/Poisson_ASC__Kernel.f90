#include "Preprocessor"

submodule ( Poisson_ASC__Form ) Poisson_ASC__Kernel

  use Basics
  use LaplacianMultipoleOld_Template
  
  implicit none

contains


  module procedure SolveCells_CSL_Kernel

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iR, &  !-- iRadius
      iE     !-- iEquation
    real ( KDR ) :: &
      R, &  !-- Radius
      S     !-- Solution
    real ( KDR ), dimension ( nAngularMomentCells ) :: &
      Zero

    Zero = 0.0_KDR

    !$OMP  parallel do &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
    !$OMP& private ( iC, iR, iE, R, S, SH_RC, SH_IC, SH_RS, SH_IS )
    do iC = 1, nCells

      if ( .not. IsProperCell ( iC ) ) &
        cycle

      call ComputeSolidHarmonicsKernel &
             ( CoordinateSystem, Center ( iC, : ), Origin, RadialEdge, &
               MaxDegree, nDimensions, GridError, SH_RC, SH_IC, SH_RS, SH_IS, &
               R, iR )

      do iE = 1, nEquations

        S = 0.0_KDR

        if ( iR > 1 .and. iR < nRadialCells ) then
          call AssembleSolutionKernel &
                 ( S, M_RC ( :, iR - 1, iE ), M_IC ( :, iR, iE ), &
                   M_RC ( :, iR, iE ), M_IC ( :, iR + 1, iE ), &
                   SH_RC, SH_IC, Delta, RadialEdge ( iR ), &
                   RadialEdge ( iR + 1 ), R, nAngularMomentCells )
        else if ( iR == 1 ) then
          call AssembleSolutionKernel &
                 ( S, Zero, M_IC ( :, iR, iE ), &
                   M_RC ( :, iR, iE ), M_IC ( :, iR + 1, iE ), &
                   SH_RC, SH_IC, Delta, RadialEdge ( iR ), &
                   RadialEdge ( iR + 1 ), R, nAngularMomentCells )
        else if ( iR == nRadialCells ) then
          call AssembleSolutionKernel &
                 ( S, M_RC ( :, iR - 1, iE ), M_IC ( :, iR, iE ), &
                   M_RC ( :, iR, iE ), Zero, &
                   SH_RC, SH_IC, Delta, RadialEdge ( iR ), &
                   RadialEdge ( iR + 1 ), R, nAngularMomentCells )
        end if !-- iR

        if ( MaxOrder > 0 ) then
          if ( iR > 1 .and. iR < nRadialCells ) then
            call AssembleSolutionKernel &
                   ( S, M_RS ( :, iR - 1, iE ), M_IS ( :, iR, iE ), &
                     M_RS ( :, iR, iE ), M_IS ( :, iR + 1, iE ), &
                     SH_RS, SH_IS, Delta, RadialEdge ( iR ), &
                     RadialEdge ( iR + 1 ), R, nAngularMomentCells )
          else if ( iR == 1 ) then
            call AssembleSolutionKernel &
                   ( S, Zero, M_IS ( :, iR, iE ), &
                     M_RS ( :, iR, iE ), M_IS ( :, iR + 1, iE ), &
                     SH_RS ( : ), SH_IS ( : ), Delta, RadialEdge ( iR ), &
                     RadialEdge ( iR + 1 ), R, nAngularMomentCells )
          else if ( iR == nRadialCells ) then
            call AssembleSolutionKernel &
                   ( S, M_RS ( :, iR - 1, iE ), M_IS ( :, iR, iE ), &
                     M_RS ( :, iR, iE ), Zero, &
                     SH_RS ( : ), SH_IS ( : ), Delta, RadialEdge ( iR ), &
                     RadialEdge ( iR + 1 ), R, nAngularMomentCells )
          end if !-- iR
        end if !-- MaxOrder

        Solution ( iC, iaSolution ( iE ) )  =  - S / FourPi

      end do !-- iE

    end do !-- iC
    !$OMP  end parallel do

  end procedure SolveCells_CSL_Kernel 


  module procedure AssembleSolutionKernel

    integer ( KDI ) :: &
      iA  !-- iAngular
    real ( KDR ) :: &
      F, &  !-- Fraction
      M_R, M_I

    F  =  ( R - R_In ) / ( R_Out - R_in )

    do iA = 1, nA

      M_R  =  F  *  M_R_In ( iA )   +   ( 1.0_KDR - F )  *  M_R_Out ( iA )
      M_I  =  F  *  M_I_In ( iA )   +   ( 1.0_KDR - F )  *  M_I_Out ( iA )

      S  =  S   +   Delta ( iA )  &
                    *  ( M_R  *  SH_I ( iA )  +  M_I  *  SH_R ( iA ) )

    end do !-- iA

  end procedure AssembleSolutionKernel


end submodule Poisson_ASC__Kernel
