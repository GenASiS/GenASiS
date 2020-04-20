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
      SolveCells_CSL_Kernel, &
      AssembleSolutionKernel

    interface

      module subroutine SolveCells_CSL_Kernel &
                          ( Solution, CoordinateSystem, IsProperCell, &
                            M_RC, M_IC, M_RS, M_IS, Center, Origin, Delta, &
                            RadialEdge, FourPi, iaSolution, MaxDegree, &
                            nDimensions, nCells, nEquations, nRadialCells, &
                            nAngularMomentCells, GridError, &
                            SH_RC, SH_IC, SH_RS, SH_IS )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          Solution
        character ( * ), intent ( in ) :: &
          CoordinateSystem
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          M_RC, M_IC, &
          M_RS, M_IS
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          Center
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Origin, &
          Delta, &
          RadialEdge
        real ( KDR ), intent ( in ) :: &
          FourPi
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSolution
        integer ( KDI ), intent ( in ) :: &
          MaxDegree, &
          nDimensions, &
          nCells, &
          nEquations, &
          nRadialCells, &
          nAngularMomentCells
        logical ( KDL ), intent ( out ) :: &
          GridError
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          SH_RC, SH_IC, &  !-- SolidHarmonics
          SH_RS, SH_IS  
      end subroutine SolveCells_CSL_Kernel

      module subroutine AssembleSolutionKernel &
                          ( S, M_R_In, M_I_In, M_R_Out, M_I_Out, SH_R, SH_I, &
                            Delta, R_In, R_Out, R, nA )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          S
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_R_In, M_I_In, &
          M_R_Out, M_I_Out, &
          SH_R, SH_I, &
          Delta
        real ( KDR ), intent ( in ) :: &
          R_In, R_Out, R
        integer ( KDI ), intent ( in ) :: &
          nA
      end subroutine AssembleSolutionKernel
      
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

    if ( associated ( Timer_SC ) ) call Timer_SC % Start ( )
     call SolveCells_CSL_Kernel &
            ( Solution % Value, C % CoordinateSystem, C % IsProperCell, &
              L % M_RC, L % M_IC, L % M_RS, L % M_IS, &
              G % Value ( :, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ), &
              L % Origin, L % Delta, L % RadialEdge, 4.0_KDR * CONSTANT % PI, &
              Solution % iaSelected, L % MaxDegree, C % nDimensions, &
              G % nValues, L % nEquations, L % nRadialCells, &
              L % nAngularMomentCells, GridError, &
              L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
              L % SolidHarmonic_RS, L % SolidHarmonic_IS )
     if ( GridError ) &
       call PROGRAM_HEADER % Abort ( )
    if ( associated ( Timer_SC ) ) call Timer_SC % Stop ( )

    if ( associated ( Timer_ES ) ) call Timer_ES % Start ( )
    call C % ExchangeGhostData ( Solution )
    if ( associated ( Timer_ES ) ) call Timer_ES % Stop ( )

    if ( associated ( Timer_BS ) ) call Timer_BS % Start ( )
    call P % Atlas % ApplyBoundaryConditionsFaces ( Solution )
    if ( associated ( Timer_BS ) ) call Timer_BS % Stop ( )

    end associate !-- L

    nullify ( G )

  end subroutine SolveMultipole_CSL


end module Poisson_ASC__Form
