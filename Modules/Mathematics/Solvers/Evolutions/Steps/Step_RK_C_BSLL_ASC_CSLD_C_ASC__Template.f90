!-- Step_RK_C_BSLL_ASC_CSLD_C_ASC is a template for a RungeKutta time step of
!   one conserved current on a bundle and another on a base space.

module Step_RK_C_BSLL_ASC_CSLD_C_ASC__Template

  !-- Step_RungeKutta_Current_BundleSingleLevelDistributed_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Operations
  use Fields
  use Step_RK_C_ASC__Template
  use Step_RK_C_ASC_1D__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_1D_Template ), abstract :: &
    Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template
      integer ( KDI ) :: &
        nFibers   = 0, &
        nSections = 0
      logical ( KDL ) :: &
        UseLimiterParameterSection
      class ( Chart_SLL_Form ), pointer :: &
        GridFiber => null ( )
      class ( Chart_SLD_Form ), pointer :: &
        GridSection => null ( )
      type ( CurrentPointerForm ), dimension ( : ), allocatable :: &
        CurrentFiber_1D, &
        CurrentSection_1D
      type ( ApplyDivergence_C_Pointer ) :: &
        ApplyDivergenceFiber, &
        ApplyDivergenceSection
      type ( ApplySources_C_Pointer ) :: &
        ApplySourcesFiber, &
        ApplySourcesSection
      type ( ApplyRelaxation_C_Pointer ) :: &
        ApplyRelaxationFiber, &
        ApplyRelaxationSection
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC
    procedure, private, pass :: &
      Compute_C_BSLL_ASC_CSLD_C_ASC
    generic, public :: &
      Compute => Compute_C_ASC_1D
    procedure, public, pass :: &
      FinalizeTemplate_C_BSLL_ASC_CSLD_C_ASC
  end type Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template

contains


  subroutine InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC &
               ( S, NameSuffix, A, B, C )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
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
      S % Type = 'a Step_RK_C_BSLL_ASC_CSLD_C_ASC' 

    call S % InitializeTemplate_C_ASC ( NameSuffix, A, B, C )

  end subroutine InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC


  subroutine Compute_C_BSLL_ASC_CSLD_C_ASC &
               ( S, Current_BSLL_ASC_CSLD, Current_ASC, Time, TimeStep, &
                 UseLimiterParameterSectionOption, UseLimiterParameterOption )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ), &
      target :: &
        Current_BSLL_ASC_CSLD
    class ( Current_ASC_Template ), intent ( inout ), target :: &
      Current_ASC
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameterSectionOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameterOption

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iS     !-- iSection

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    S % UseLimiterParameter = .true.
    if ( present ( UseLimiterParameterOption ) ) &
      S % UseLimiterParameter = UseLimiterParameterOption

    S % UseLimiterParameterSection = .true.
    if ( present ( UseLimiterParameterSectionOption ) ) &
      S % UseLimiterParameterSection = UseLimiterParameterSectionOption

    select type ( Chart => Current_ASC % Atlas_SC % Chart )
    class is ( Chart_SL_Template )
      S % Grid => Chart
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_C_ASC', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart
    S % Current => Current_ASC % Current ( )

    associate &
      ( CB => Current_BSLL_ASC_CSLD, &
        B  => Current_BSLL_ASC_CSLD % Bundle_SLL_ASC_CSLD )

    S % nFibers   =  B % nFibers
    S % GridFiber => B % Fiber_CSLL
    allocate ( S % CurrentFiber_1D ( S % nFibers ) )
    do iF = 1, S % nFibers
      select type ( CA => CB % Fiber % Atlas ( iF ) % Element )
      class is ( Current_ASC_Template )
        S % CurrentFiber_1D ( iF ) % Pointer => CA % Current ( )
      end select !-- CA
    end do !-- iF

    S % nSections   =  CB % nSections
    S % GridSection => B  % Base_CSLD
    allocate ( S % CurrentSection_1D ( S % nSections ) )
    do iS = 1, S % nSections
      associate &
        ( CA => CB % Section_ASC ( iS ) % Element )
      S % CurrentSection_1D ( iS ) % Pointer => CA % Current ( )
      end associate !-- CA
    end do !-- iS

    end associate !-- CurrentBundle, etc.

    ! call AllocateStorage ( S )
    ! call S % LoadSolution ( S % Solution_1D, S % Current_1D )

    ! call S % ComputeTemplate ( Time, TimeStep )

    ! call S % StoreSolution ( S % Current_1D, S % Solution_1D )
    ! call DeallocateStorage ( S )

    deallocate ( S % CurrentSection_1D )
    nullify ( S % GridSection )
    S % nSections = 0

    deallocate ( S % CurrentSection_1D )
    nullify ( S % GridSection )
    S % nFibers =  0

    S % Current => null ( )
    S % Grid    => null ( )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_C_BSLL_ASC_CSLD_C_ASC


  impure elemental subroutine FinalizeTemplate_C_BSLL_ASC_CSLD_C_ASC ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    call S % FinalizeTemplate_C_ASC ( )

  end subroutine FinalizeTemplate_C_BSLL_ASC_CSLD_C_ASC


end module Step_RK_C_BSLL_ASC_CSLD_C_ASC__Template
