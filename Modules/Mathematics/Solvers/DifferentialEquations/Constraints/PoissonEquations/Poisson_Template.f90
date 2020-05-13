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
      iTimerAssembleSolution = 0, &
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
    procedure ( S ), public, pass, deferred :: &
      Solve
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, public, pass :: &
      AssembleSolution
    procedure ( ASC ), private, pass, deferred :: &
      AssembleSolutionContributions
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

    subroutine ASC ( P, Solution, Delta_M_FourPi, iA, iSH_0 )
      use Basics
      import PoissonTemplate
      implicit none
      class ( PoissonTemplate ), intent ( inout ) :: &
        P
      class ( * ), intent ( inout ) :: &
        Solution
      real ( KDR ), intent ( in ) :: &
        Delta_M_FourPi
      integer ( KDI ), intent ( in ) :: &
        iA, &  
        iSH_0  
    end subroutine ASC

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
           ( 'AssembleSolution', P % iTimerAssembleSolution, &
             Level = BaseLevel + 1 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ExchangeSolution', P % iTimerExchangeSolution, &
             Level = BaseLevel + 1 )
    call PROGRAM_HEADER % AddTimer &
           ( 'BoundarySolution', P % iTimerBoundarySolution, &
             Level = BaseLevel + 1 )

  end subroutine InitializeTimers


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


  subroutine AssembleSolution ( P, Solution )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( * ), intent ( inout ) :: &
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
    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( P % LaplacianMultipole) ) then
      call Show ( 'LaplacianMultipole not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_Template', 'module', CONSOLE % ERROR )
      call Show ( 'AssembleSolution', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    Timer  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerAssembleSolution )
    if ( associated ( Timer ) ) call Timer % Start ( )

    associate ( L => P % LaplacianMultipole )

    iSH_0  =>  iSH ( 1 )
    iSH_1  =>  iSH ( 2 )
    iSH_2  =>  iSH ( 3 )

        iA  =  1
       iSH  =  [ 1, 2, 3 ]
    iSH_PD  =  4

    call L % Clear ( Solution )

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

        call P % AssembleSolutionContributions &
               ( Solution, Delta_M_FourPi, iA, iSH_0 )

         iA  =  iA + 1
        iSH  =  cshift ( iSH, -1 )

      end do !-- iL
    end do !-- iM

    end associate !-- L

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine AssembleSolution


end module Poisson_Template
