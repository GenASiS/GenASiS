module Poisson_Template

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, abstract :: PoissonTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerSolve = 0, &
      iTimerSolveCells = 0, &
      iTimerExchangeSolution = 0, &
      iTimerBoundarySolution = 0, &
      nEquations = 0, &
      MaxDegree = 0
    character ( LDF ) :: &
      Type = '', &
      Name = '', &
      SolverType = ''
    class ( LaplacianMultipoleTemplate ), allocatable :: &
      LaplacianMultipole
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      InitializeTimers
    procedure ( S ), public, pass, deferred :: &
      Solve
    procedure, public, pass :: &
      FinalizeTemplate
  end type PoissonTemplate

  abstract interface 

    subroutine S ( P, Solution, Source )
      use Basics
      use Manifolds
      import PoissonTemplate
      class ( PoissonTemplate ), intent ( inout ) :: &
        P
      class ( FieldAtlasTemplate ), intent ( inout ) :: &
        Solution
      class ( FieldAtlasTemplate ), intent ( in ) :: &
        Source
    end subroutine S

  end interface

contains


  subroutine InitializeTemplate &
               ( P, A, SolverType, MaxDegreeOption, nEquationsOption )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      SolverType
    integer ( KDI ), intent ( in ), optional :: &
      MaxDegreeOption, &
      nEquationsOption

    P % IGNORABILITY = A % IGNORABILITY

    if ( P % Type == '' ) &
      P % Type = 'a Poisson' 

    P % Name = 'Poisson_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )

    P % nEquations = 1
    if ( present ( nEquationsOption ) ) &
      P % nEquations = nEquationsOption

    P % SolverType = SolverType

    P % MaxDegree = 0
    if ( present ( MaxDegreeOption ) ) &
      P % MaxDegree = MaxDegreeOption

  end subroutine InitializeTemplate


  subroutine InitializeTimers ( P, BaseLevel )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    call PROGRAM_HEADER % AddTimer &
           ( 'PoissonSolve', P % iTimerSolve, Level = BaseLevel )

    if ( allocated ( P % LaplacianMultipole ) ) then

      associate ( LM => P % LaplacianMultipole )
      call LM % InitializeTimers ( BaseLevel + 1 )
      end associate !-- LM

      call PROGRAM_HEADER % AddTimer &
             ( 'SolveCells', P % iTimerSolveCells, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ExchangeSolution', P % iTimerExchangeSolution, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'BoundarySolution', P % iTimerBoundarySolution, &
               Level = BaseLevel + 1 )

    end if
  end subroutine InitializeTimers


  impure elemental subroutine FinalizeTemplate ( P )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P

    if ( allocated ( P % LaplacianMultipole ) ) &
      deallocate ( P % LaplacianMultipole )

    if ( P % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    
  end subroutine FinalizeTemplate


end module Poisson_Template
