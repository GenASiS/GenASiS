module Poisson_H__Form

  !-- Poisson_Header__Form

  use Basics
  use Fields
  use Laplacian_M_H__Form
  
  implicit none
  private

  type, public :: Poisson_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
!       iTimerSolve = 0, &
!       iTimerCombineMoments = 0, &
!       iTimerClearSolution = 0, &
!       iTimerLocalSolution = 0, &
!       iTimerExchangeSolution = 0, &
!       iTimerBoundarySolution = 0, &
      nEquations = 0, &
      MaxDegree = 0
    character ( LDF ) :: &
      Type = '', &
      Name = '', &
      SolverType = ''
    class ( Laplacian_M_H_Form ), allocatable :: &
      Laplacian_M
  contains
    procedure, public, pass :: &
      Initialize_H
!     procedure, public, pass :: &
!       InitializeTimers
    procedure, public, pass :: &
      Show => Show_P
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
    procedure, private, pass :: &
      Solve_M
    procedure, private, pass :: &
      CombineMoments
    procedure, private, pass :: &
      CombineMomentsLocal
    procedure, private, pass :: &
      ExchangeSolution
    procedure, private, pass :: &
      ApplyBoundarySolution
  end type Poisson_H_Form


contains


  subroutine Initialize_H &
               ( P, GA, SolverType, MaxDegreeOption, nEquationsOption )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( Geometry_F_A_Form ), intent ( in ) :: &
      GA
    character ( * ), intent ( in ) :: &
      SolverType
    integer ( KDI ), intent ( in ), optional :: &
      MaxDegreeOption, &
      nEquationsOption

    P % IGNORABILITY  =  GA % IGNORABILITY

    if ( P % Type == '' ) &
      P % Type = 'a Poisson' 

    P % Name = 'Poisson'

    call Show ( 'Initializing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )

    P % nEquations  =  1
    if ( present ( nEquationsOption ) ) &
      P % nEquations = nEquationsOption

    P % SolverType  =  SolverType

    P % MaxDegree  =  0
    if ( present ( MaxDegreeOption ) ) &
      P % MaxDegree  =  MaxDegreeOption

  end subroutine Initialize_H


!   subroutine InitializeTimers ( P, BaseLevel )

!     class ( PoissonTemplate ), intent ( inout ) :: &
!       P
!     integer ( KDI ), intent ( in ) :: &
!       BaseLevel

!     call PROGRAM_HEADER % AddTimer &
!            ( 'PoissonSolve', P % iTimerSolve, Level = BaseLevel )

!     if ( allocated ( P % LaplacianMultipoleOld_1 ) ) then
!       associate ( L => P % LaplacianMultipoleOld_1 )
!       call L % InitializeTimers ( BaseLevel + 1 )
!       end associate !-- L
!     end if

!     if ( allocated ( P % LaplacianMultipoleOld_2 ) ) then
!       associate ( L => P % LaplacianMultipoleOld_2 )
!       call L % InitializeTimers ( BaseLevel + 1 )
!       end associate !-- L
!     end if

!     if ( allocated ( P % LaplacianMultipole ) ) then
!       associate ( L => P % LaplacianMultipole )
!       call L % InitializeTimers ( BaseLevel + 1 )
!       end associate !-- L
!     end if

!     call PROGRAM_HEADER % AddTimer &
!            ( 'CombineMoments', P % iTimerCombineMoments, &
!              Level = BaseLevel + 1 )
!       call PROGRAM_HEADER % AddTimer &
!              ( 'ClearSolution', P % iTimerClearSolution, &
!                Level = BaseLevel + 2 )
!       call PROGRAM_HEADER % AddTimer &
!              ( 'LocalSolution', P % iTimerLocalSolution, &
!                Level = BaseLevel + 2 )
!       call PROGRAM_HEADER % AddTimer &
!              ( 'ExchangeSolution', P % iTimerExchangeSolution, &
!                Level = BaseLevel + 2 )
!       call PROGRAM_HEADER % AddTimer &
!              ( 'BoundarySolution', P % iTimerBoundarySolution, &
!                Level = BaseLevel + 2 )

!   end subroutine InitializeTimers


  subroutine Show_P ( P )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( P % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    call Show ( P % SolverType, 'SolverType', P % IGNORABILITY )

    if ( allocated ( P % Laplacian_M ) ) &
      call P % Laplacian_M % Show ( )

  end subroutine Show_P


  subroutine Solve ( P, Solution_A, Source_A )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Solution_A
    class ( FieldSet_A_Form ), intent ( in ) :: &
      Source_A

!     type ( TimerForm ), pointer :: &
!       Timer

!     Timer  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerSolve )
!     if ( associated ( Timer ) ) call Timer % Start ( )

    select case ( trim ( P % SolverType ) )
    case ( 'MULTIPOLE' )

      call P % Solve_M ( Solution_A, Source_A )

    case default
      call Show ( 'Solver type not supported', CONSOLE % ERROR )
      call Show ( P % SolverType, 'Type', CONSOLE % ERROR )
      call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SolverType

!     if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Solve


  impure elemental subroutine Finalize ( P )

    type ( Poisson_H_Form ), intent ( inout ) :: &
      P

    if ( allocated ( P % Laplacian_M ) ) &
      deallocate ( P % Laplacian_M )

    if ( P % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    
  end subroutine Finalize


  subroutine Solve_M ( P, Solution_A, Source_A )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Solution_A
    class ( FieldSet_A_Form ), intent ( in ) :: &
      Source_A

    call Show ( 'Poisson solve, multipole', P % IGNORABILITY + 2 )
    call Show ( P % Name, 'Name', P % IGNORABILITY + 2 )

    if ( allocated ( P % Laplacian_M ) ) then
      associate ( L  =>  P % Laplacian_M )
      call L % ComputeMoments ( Source_A )
      call P % CombineMoments ( Solution_A )
      end associate !-- LA
    else
      call Show ( 'Laplacian_M not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Solve_M', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine Solve_M


  subroutine CombineMoments ( P, Solution_A )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Solution_A

!     type ( TimerForm ), pointer :: &
!       Timer, &
!       Timer_CS, &
!       Timer_LS, &
!       Timer_ES, &
!       Timer_BS

    if ( .not. allocated ( P % Laplacian_M ) ) then
      call Show ( 'Laplacian_M not allocated', CONSOLE % ERROR )
      call Show ( 'Poisson_H_Form', 'module', CONSOLE % ERROR )
      call Show ( 'CombineMoments', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

!     Timer     =>  PROGRAM_HEADER % TimerPointer ( P % iTimerCombineMoments )
!     Timer_CS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerClearSolution )
!     Timer_LS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerLocalSolution )
!     Timer_ES  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerExchangeSolution )
!     Timer_BS  =>  PROGRAM_HEADER % TimerPointer ( P % iTimerBoundarySolution )

!     if ( associated ( Timer ) ) call Timer % Start ( )

    call Show ( 'Combining Moments', P % IGNORABILITY + 2 )

!     if ( associated ( Timer_CS ) ) call Timer_CS % Start ( )
    call Solution_A % Clear ( )
!     if ( associated ( Timer_CS ) ) call Timer_CS % Stop ( )

!     if ( associated ( Timer_LS ) ) call Timer_LS % Start ( )
!     if ( allocated ( P % LaplacianMultipoleOld_2 ) ) then
!       call P % CombineMomentsLocalOld_2 ( Solution )
    if ( allocated ( P % Laplacian_M ) ) then
      call P % CombineMomentsLocal ( Solution_A )
    end if
!     if ( associated ( Timer_LS ) ) call Timer_LS % Stop ( )

!     if ( associated ( Timer_ES ) ) call Timer_ES % Start ( )
    call P % ExchangeSolution ( Solution_A )
!     if ( associated ( Timer_ES ) ) call Timer_ES % Stop ( )

!     if ( associated ( Timer_BS ) ) call Timer_BS % Start ( )
    call P % ApplyBoundarySolution ( Solution_A )
!     if ( associated ( Timer_BS ) ) call Timer_BS % Stop ( )

!     if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine CombineMoments


  subroutine CombineMomentsLocal ( P, Solution_A )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Solution_A

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'CombineMomentsLocal', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine CombineMomentsLocal


  subroutine ExchangeSolution ( P, Solution_A )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Solution_A

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ExchangeSolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ExchangeSolution


  subroutine ApplyBoundarySolution ( P, Solution_A )

    class ( Poisson_H_Form ), intent ( inout ) :: &
      P
    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Solution_A

    call Show ( 'Subroutine should be overidden', CONSOLE % ERROR )
    call Show ( 'Poisson_H__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ApplyBoundarySolution', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ApplyBoundarySolution


end module Poisson_H__Form
