module Poisson_ASC__Form

  !-- Poisson_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template
  use LaplacianMultipole_ASC__Form
  use Poisson_Template

  implicit none
  private

  type, public, extends ( PoissonTemplate ) :: Poisson_ASC_Form
    class ( Atlas_SC_Form ), pointer :: &
      Atlas => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
  end type Poisson_ASC_Form

    private :: &
      SolveMultipole_CSL

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      SolveCells_CSL_Kernel

    interface

      module subroutine SolveCells_CSL_Kernel &
                          ( CoordinateSystem, IsProperCell, Center, &
                            nDimensions, nCells, ComputeSolidHarmonics )
        use Basics
        character ( * ), intent ( in ) :: &
          CoordinateSystem
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          Center
        integer ( KDI ), intent ( in ) :: &
          nDimensions, &
          nCells
        procedure ( ), pointer :: &
          ComputeSolidHarmonics
      end subroutine SolveCells_CSL_Kernel

    end interface

contains


  subroutine Initialize &
               ( P, A, SolverType, MaxDegreeOption, nEquationsOption )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      SolverType
    integer ( KDI ), intent ( in ), optional :: &
      MaxDegreeOption, &
      nEquationsOption

    if ( P % Type == '' ) &
      P % Type = 'a Poisson_ASC' 

    call P % InitializeTemplate &
           ( A, SolverType, MaxDegreeOption, nEquationsOption )

    P % Atlas => A

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE' )
      allocate ( LaplacianMultipole_ASC_Form :: P % LaplacianMultipole )
      select type ( L => P % LaplacianMultipole )
      class is ( LaplacianMultipole_ASC_Form )
        call L % Initialize ( A, P % MaxDegree, P % nEquations )
      end select !-- L
    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine Initialize


  subroutine Solve ( P, Solution, Source )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source

    class ( StorageForm ), pointer :: &
      Source_S, &
      Solution_S
    type ( TimerForm ), pointer :: &
      Timer

    Timer  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerSolve )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( Source )
    class is ( Storage_ASC_Form )
    Source_S => Source % Storage ( )

    select type ( Solution )
    class is ( Storage_ASC_Form )
    Solution_S => Solution % Storage ( )

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE' )
      select type ( C => P % Atlas % Chart )
      class is ( Chart_SLD_Form )
        call SolveMultipole_CSL ( P, C, Solution_S, Source_S )
      class default
        call Show ( 'Chart type not supported', CONSOLE % ERROR )
        call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
        call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( P % SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SolverType

    end select !-- Solution
    end select !-- Source

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Solve


  impure elemental subroutine Finalize ( P )

    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P

    call P % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SolveMultipole_CSL ( P, C, Solution, Source )
 
    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( StorageForm ), intent ( inout ) :: &
      Solution
    class ( StorageForm ), intent ( in ) :: &
      Source

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iR, &  !-- iRadius
      iE     !-- iEquation
    real ( KDR ) :: &
      R, &  !-- Radius
      S     !-- Solution
    real ( KDR ), dimension ( : ), allocatable :: &
      Zero
    logical ( KDL ) :: &
      GridError
    type ( TimerForm ), pointer :: &
      Timer_SC, &
      Timer_ES, &
      Timer_BS
    class ( GeometryFlatForm ), pointer :: &
      G

    Timer_SC  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerSolveCells )
    Timer_ES  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerExchangeSolution )
    Timer_BS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerBoundarySolution )

    call Show ( 'Poisson solve, multipole', P % IGNORABILITY + 2 )
    call Show ( P % Name, 'Name', P % IGNORABILITY + 2 )

    associate ( L  =>  P % LaplacianMultipole )

    call L % ComputeMoments ( Source )

    G => C % Geometry ( )

    allocate ( Zero ( L % nAngularMomentCells ) )
    Zero = 0.0_KDR

    do iE = 1, L % nEquations

      if ( associated ( Timer_SC ) ) call Timer_SC % Start ( )

!      call SolveCells_CSL_Kernel &
!             ( C % CoordinateSystem, C % IsProperCell, &
!               G % Value ( :, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ), &
!               C % nDimensions, G % nValues, L % ComputeSolidHarmonicsKernel )
!      if ( GridError ) &
!        call PROGRAM_HEADER % Abort ( )

!      !$OMP parallel do private ( iC, iR, R, S )
      do iC = 1, G % nValues

        if ( .not. C % IsProperCell ( iC ) ) &
          cycle

        call ComputeSolidHarmonicsKernel &
               ( C % CoordinateSystem, &
                 G % Value ( iC, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ), &
                 L % Origin, L % RadialEdge, L % MaxDegree, C % nDimensions, &
                 GridError, L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
                 L % SolidHarmonic_RS, L % SolidHarmonic_IS, R, iR )
        if ( GridError ) &
          call PROGRAM_HEADER % Abort ( )

        S = 0.0_KDR

        if ( iR > 1 .and. iR < L % nRadialCells ) then
          call SolveKernel &
                 ( S, &
                   L % M_RC ( :, iR - 1, iE ), L % M_IC ( :, iR, iE ), &
                   L % M_RC ( :, iR, iE ), L % M_IC ( :, iR + 1, iE ), &
                   L % SolidHarmonic_RC ( : ), L % SolidHarmonic_IC ( : ), &
                   L % Delta, L % RadialEdge ( iR ), &
                   L % RadialEdge ( iR + 1 ), R )
        else if ( iR == 1 ) then
          call SolveKernel &
                 ( S, &
                   Zero, L % M_IC ( :, iR, iE ), &
                   L % M_RC ( :, iR, iE ), L % M_IC ( :, iR + 1, iE ), &
                   L % SolidHarmonic_RC ( : ), L % SolidHarmonic_IC ( : ), &
                   L % Delta, L % RadialEdge ( iR ), &
                   L % RadialEdge ( iR + 1 ), R )
        else if ( iR == L % nRadialCells ) then
          call SolveKernel &
                 ( S, &
                   L % M_RC ( :, iR - 1, iE ), L % M_IC ( :, iR, iE ), &
                   L % M_RC ( :, iR, iE ), Zero, &
                   L % SolidHarmonic_RC ( : ), L % SolidHarmonic_IC ( : ), &
                   L % Delta, L % RadialEdge ( iR ), &
                   L % RadialEdge ( iR + 1 ), R )
        end if !-- iR

        if ( L % MaxOrder > 0 ) then
          if ( iR > 1 .and. iR < L % nRadialCells ) then
            call SolveKernel &
                   ( S, &
                     L % M_RS ( :, iR - 1, iE ), L % M_IS ( :, iR, iE ), &
                     L % M_RS ( :, iR, iE ), L % M_IS ( :, iR + 1, iE ), &
                     L % SolidHarmonic_RS ( : ), L % SolidHarmonic_IS ( : ), &
                     L % Delta, L % RadialEdge ( iR ), &
                     L % RadialEdge ( iR + 1 ), R )
          else if ( iR == 1 ) then
            call SolveKernel &
                   ( S, &
                     Zero, L % M_IS ( :, iR, iE ), &
                     L % M_RS ( :, iR, iE ), L % M_IS ( :, iR + 1, iE ), &
                     L % SolidHarmonic_RS ( : ), L % SolidHarmonic_IS ( : ), &
                     L % Delta, L % RadialEdge ( iR ), &
                     L % RadialEdge ( iR + 1 ), R )
          else if ( iR == L % nRadialCells ) then
            call SolveKernel &
                   ( S, &
                     L % M_RS ( :, iR - 1, iE ), L % M_IS ( :, iR, iE ), &
                     L % M_RS ( :, iR, iE ), Zero, &
                     L % SolidHarmonic_RS ( : ), L % SolidHarmonic_IS ( : ), &
                     L % Delta, L % RadialEdge ( iR ), &
                     L % RadialEdge ( iR + 1 ), R )
          end if !-- iR
        end if !-- MaxOrder

        Solution % Value ( iC, Solution % iaSelected ( iE ) )  &
          =  - S / ( 4.0_KDR  *  CONSTANT % PI )

      end do !-- iC
!      !$OMP end parallel do

      if ( associated ( Timer_SC ) ) call Timer_SC % Stop ( )

    end do !-- iE

    if ( associated ( Timer_ES ) ) call Timer_ES % Start ( )
    call C % ExchangeGhostData ( Solution )
    if ( associated ( Timer_ES ) ) call Timer_ES % Stop ( )

    if ( associated ( Timer_BS ) ) call Timer_BS % Start ( )
    call P % Atlas % ApplyBoundaryConditionsFaces ( Solution )
    if ( associated ( Timer_BS ) ) call Timer_BS % Stop ( )

    end associate !-- L

    nullify ( G )

  end subroutine SolveMultipole_CSL


  subroutine SolveKernel &
               ( S, M_R_In, M_I_In, M_R_Out, M_I_Out, SH_R, SH_I, Delta, &
                 R_In, R_Out, R )

    real ( KDR ), intent ( inout ) :: &
      S
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_R_In, M_I_In, &
      M_R_Out, M_I_Out, &
      SH_R, SH_I, &
      Delta
    real ( KDR ), intent ( in ) :: &
      R_In, R_Out, R

    integer ( KDI ) :: &
      iA  !-- iAngular
    real ( KDR ) :: &
      F, &  !-- Fraction
      M_R, M_I

    F  =  ( R - R_In ) / ( R_Out - R_in )

    do iA = 1, size ( Delta )

      M_R  =  F  *  M_R_In ( iA )   +   ( 1.0_KDR - F )  *  M_R_Out ( iA )
      M_I  =  F  *  M_I_In ( iA )   +   ( 1.0_KDR - F )  *  M_I_Out ( iA )

      S  =  S   +   Delta ( iA )  &
                    *  ( M_R  *  SH_I ( iA )  +  M_I  *  SH_R ( iA ) )

    end do

  end subroutine SolveKernel


end module Poisson_ASC__Form
