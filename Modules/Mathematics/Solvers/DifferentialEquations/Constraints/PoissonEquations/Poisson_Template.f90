module Poisson_Template

  use Basics
  use Manifolds
  use LaplacianMultipoleOld_Template
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, abstract :: PoissonTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerSolve = 0, &
      iTimerCombineMoments = 0, &
      iTimerClearSolution = 0, &
      iTimerLocalSolution = 0, &
      iTimerExchangeSolution = 0, &
      iTimerBoundarySolution = 0, &
      nEquations = 0, &
      MaxDegree = 0
    character ( LDF ) :: &
      Type = '', &
      Name = '', &
      SolverType = ''
    class ( LaplacianMultipoleOldTemplate ), allocatable :: &
      LaplacianMultipoleOld
    class ( LaplacianMultipoleTemplate ), allocatable :: &
      LaplacianMultipole
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      InitializeTimers
    procedure, public, pass :: &
      Solve
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( SO ), private, pass, deferred :: &
      SolveOld
    procedure, private, pass :: &
      SolveMultipole
    procedure, private, pass :: &
      CombineMoments
    procedure, private, pass :: &
      CombineMomentsLocal
    procedure ( CMA ), private, pass, deferred :: &
      CombineMomentAtlas
    procedure ( ES ), private, pass, deferred :: &
      ExchangeSolution
    procedure ( ABS ), private, pass, deferred :: &
      ApplyBoundarySolution
  end type PoissonTemplate

  abstract interface 

    subroutine SO ( P, Solution, Source )
      use Manifolds
      import PoissonTemplate
      implicit none
      class ( PoissonTemplate ), intent ( inout ) :: &
        P
      class ( FieldAtlasTemplate ), intent ( inout ) :: &
        Solution
      class ( FieldAtlasTemplate ), intent ( in ) :: &
        Source
    end subroutine SO

    subroutine CMA ( P, Solution, Delta_M_FourPi, iA, iSH_0 )
      use Basics
      use Manifolds
      import PoissonTemplate
      implicit none
      class ( PoissonTemplate ), intent ( inout ) :: &
        P
      class ( FieldAtlasTemplate ), intent ( inout ) :: &
        Solution
      real ( KDR ), intent ( in ) :: &
        Delta_M_FourPi
      integer ( KDI ), intent ( in ) :: &
        iA, &  
        iSH_0  
    end subroutine CMA

    subroutine ES ( P, Solution )
      use Manifolds
      import PoissonTemplate
      implicit none
      class ( PoissonTemplate ), intent ( inout ) :: &
        P
      class ( FieldAtlasTemplate ), intent ( inout ) :: &
        Solution
    end subroutine ES

    subroutine ABS ( P, Solution )
      use Manifolds
      import PoissonTemplate
      implicit none
      class ( PoissonTemplate ), intent ( inout ) :: &
        P
      class ( FieldAtlasTemplate ), intent ( inout ) :: &
        Solution
    end subroutine ABS

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

    if ( allocated ( P % LaplacianMultipoleOld ) ) then
      associate ( L => P % LaplacianMultipoleOld )
      call L % InitializeTimers ( BaseLevel + 1 )
      end associate !-- L
    end if

    if ( allocated ( P % LaplacianMultipole ) ) then
      associate ( L => P % LaplacianMultipole )
      call L % InitializeTimers ( BaseLevel + 1 )
      end associate !-- L
    end if

    call PROGRAM_HEADER % AddTimer &
           ( 'CombineMoments', P % iTimerCombineMoments, &
             Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ClearSolution', P % iTimerClearSolution, &
               Level = BaseLevel + 2 )
      call PROGRAM_HEADER % AddTimer &
             ( 'LocalSolution', P % iTimerLocalSolution, &
               Level = BaseLevel + 2 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ExchangeSolution', P % iTimerExchangeSolution, &
               Level = BaseLevel + 2 )
      call PROGRAM_HEADER % AddTimer &
             ( 'BoundarySolution', P % iTimerBoundarySolution, &
               Level = BaseLevel + 2 )

  end subroutine InitializeTimers


  subroutine Solve ( P, Solution, Source )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source

    type ( TimerForm ), pointer :: &
      Timer

    Timer  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerSolve )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE_OLD', 'MULTIPOLE' )

      call P % SolveMultipole ( Solution, Source )

    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( P % SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Poisson_Template', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SolverType

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Solve


  impure elemental subroutine FinalizeTemplate ( P )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P

    if ( allocated ( P % LaplacianMultipole ) ) &
      deallocate ( P % LaplacianMultipole )
    if ( allocated ( P % LaplacianMultipoleOld ) ) &
      deallocate ( P % LaplacianMultipoleOld )

    if ( P % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    
  end subroutine FinalizeTemplate


  subroutine SolveMultipole ( P, Solution, Source )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source

    call Show ( 'Poisson solve, multipole', P % IGNORABILITY + 2 )
    call Show ( P % Name, 'Name', P % IGNORABILITY + 2 )

    if ( allocated ( P % LaplacianMultipoleOld ) ) then
      call P % SolveOld ( Solution, Source )
    else if ( allocated ( P % LaplacianMultipole ) ) then
      associate ( L  =>  P % LaplacianMultipole )
        call L % ComputeMoments ( Source )
        call P % CombineMoments ( Solution )
      end associate !-- L
    else
      call Show ( 'LaplacianMultipole not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_Template', 'module', CONSOLE % ERROR )
      call Show ( 'SolveMultipole', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine SolveMultipole


  subroutine CombineMoments ( P, Solution )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution

    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_CS, &
      Timer_LS, &
      Timer_ES, &
      Timer_BS

    if ( .not. allocated ( P % LaplacianMultipole) ) then
      call Show ( 'LaplacianMultipole not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_Template', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMoments', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    Timer     =>  PROGRAM_HEADER % TimerPointer ( P % iTimerCombineMoments )
    Timer_CS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerClearSolution )
    Timer_LS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerLocalSolution )
    Timer_ES  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerExchangeSolution )
    Timer_BS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerBoundarySolution )

    if ( associated ( Timer ) ) call Timer % Start ( )

    if ( associated ( Timer_CS ) ) call Timer_CS % Start ( )
    call Solution % Clear ( )
    if ( associated ( Timer_CS ) ) call Timer_CS % Stop ( )

    if ( associated ( Timer_LS ) ) call Timer_LS % Start ( )
    call P % CombineMomentsLocal ( Solution )
    if ( associated ( Timer_LS ) ) call Timer_LS % Stop ( )

    if ( associated ( Timer_ES ) ) call Timer_ES % Start ( )
    call P % ExchangeSolution ( Solution )
    if ( associated ( Timer_ES ) ) call Timer_ES % Stop ( )

    if ( associated ( Timer_BS ) ) call Timer_BS % Start ( )
    call P % ApplyBoundarySolution ( Solution )
    if ( associated ( Timer_BS ) ) call Timer_BS % Stop ( )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine CombineMoments


  subroutine CombineMomentsLocal ( P, Solution )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      Solution

    integer ( KDI ) :: &
      iA, &   !-- iAngularMoment
      iM, &   !-- iOrder
      iL, &   !-- iDegree
      iSH_PD  !-- iSolidHarmonic_PreviousDiagonal
    integer ( KDI ), pointer :: &
      iSH_0, &  !-- iSolidHarmonic_Current
      iSH_1, &  !-- iSolidHarmonic_Previous_1
      iSH_2     !-- iSolidHarmonic_Previous_2
    integer ( KDI ), dimension ( 3 ), target :: &
      iSH
    real ( KDR ) :: &
      Delta_M_FourPi

    associate ( L => P % LaplacianMultipole )

    iSH_0  =>  iSH ( 1 )
    iSH_1  =>  iSH ( 2 )
    iSH_2  =>  iSH ( 3 )

        iA  =  1
       iSH  =  [ 1, 2, 3 ]
    iSH_PD  =  4

    do iM  =  0, L % MaxOrder

      !-- ( L, M ) = ( iM, iM )
      !-- Note iL = iM
      if ( iM  ==  0 ) then
        Delta_M_FourPi  =  - 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI )
        call L % ComputeSolidHarmonics_0_0 ( iSH_0, iSH_PD )
      else
        Delta_M_FourPi  =  - 2.0_KDR / ( 4.0_KDR  *  CONSTANT % PI )
        call L % ComputeSolidHarmonics_iM_iM ( iM, iSH_0, iSH_PD )
      end if

      do iL  =  iM, L % MaxDegree

        if ( iL  ==  iM + 1 ) then
          !-- ( L, M ) = ( iM + 1, iM )
          !-- Note iL = iM + 1
          call L % ComputeSolidHarmonics_iL_iM_1 &
                 ( iM, iSH_0, iSH_1 )
        else if ( iL  >=  iM  +  2 ) then
          !-- ( L, M ) = ( iL, iM )
          call L % ComputeSolidHarmonics_iL_iM_2 &
                 ( iL, iM, iSH_0, iSH_1, iSH_2 )
        end if

        call P % CombineMomentAtlas &
               ( Solution, Delta_M_FourPi, iA, iSH_0 )

         iA  =  iA + 1
        iSH  =  cshift ( iSH, -1 )

      end do !-- iL
    end do !-- iM

    end associate !-- L

  end subroutine CombineMomentsLocal


end module Poisson_Template
