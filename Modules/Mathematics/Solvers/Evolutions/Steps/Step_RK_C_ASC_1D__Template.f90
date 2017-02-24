!-- Step_RK_C_ASC is a template for a RungeKutta time step of an array of
!   conserved currents.

module Step_RK_C_ASC_1D__Template

  !-- Step_RungeKutta_Current_AtlasSingleChart_1D_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Fields
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ), abstract :: &
    Step_RK_C_ASC_1D_Template
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_1D
    procedure, private, pass :: &
      Compute_C_1D
    generic, public :: &
      Compute => Compute_C_1D
    procedure, public, pass :: &
      FinalizeTemplate_C_1D
  end type Step_RK_C_ASC_1D_Template

contains


  subroutine InitializeTemplate_C_1D ( S, NameSuffix, A, B, C )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_ASC_1D' 

    call S % InitializeTemplate_C ( NameSuffix, A, B, C )

  end subroutine InitializeTemplate_C_1D


  subroutine Compute_C_1D &
               ( S, Current_ASC_1D, Grid, Time, TimeStep, &
                 UseLimiterParameterOption )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    class ( Current_ASC_ElementForm ), dimension ( : ), intent ( inout ), &
      target :: &
        Current_ASC_1D
    class ( * ), intent ( in ), target :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    logical ( KDL ), dimension ( : ), intent ( in ), optional :: &
      UseLimiterParameterOption

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    ! S % UseLimiterParameter_C = .true.
    ! if ( present ( UseLimiterParameterOption ) ) &
    !   S % UseLimiterParameter_C = UseLimiterParameterOption

    ! S % Grid      => Grid
    ! S % Current_C => Current_ASC % Current ( )
    ! call AllocateStorage ( S, S % Current_C )
    ! call S % LoadSolution ( S % Solution_C, S % Current_C )

    ! call S % ComputeTemplate ( Time, TimeStep )

    ! call S % StoreSolution ( S % Current_C, S % Solution_C )
    ! call DeallocateStorage ( S )
    ! S % Current_C => null ( )
    ! S % Grid      => null ( )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_C_1D


  impure elemental subroutine FinalizeTemplate_C_1D ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    call S % FinalizeTemplate_C ( )

  end subroutine FinalizeTemplate_C_1D


end module Step_RK_C_ASC_1D__Template
