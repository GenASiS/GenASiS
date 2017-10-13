!-- Step_RK_C_BSLL_ASC_CSLD is a template for a RungeKutta time step of
!   a conserved current on a bundle.

module Step_RK_C_BSLL_ASC_CSLD__Template

  !-- Step_RungeKutta_Current_BundleSingleLevelDistributed_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Operations
  use Fields
  use Increments
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ), abstract :: &
    Step_RK_C_BSLL_ASC_CSLD_Template
      integer ( KDI ) :: &
        iTimerSections = 0, &
        iTimerStoreSections = 0, &
        iTimerFibers = 0, &
        iTimerLoadSections = 0, &
        nFibers   = 0, &
        nSections = 0
      type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
        BoundaryFluence_CSL_S
      type ( VariableGroupForm ), dimension ( : ), allocatable :: &
        Solution_BSLL_ASC_CSLD_S, &
        Y_BSLL_ASC_CSLD_S
      type ( Storage_BSLL_ASC_CSLD_Form ), dimension ( : ), allocatable :: &
        K_BSLL_ASC_CSLD
      class ( Chart_SLL_Form ), pointer :: &
        Chart_F => null ( )
      class ( Bundle_SLL_ASC_CSLD_Form ), pointer :: &
        Bundle_SLL_ASC_CSLD => null ( )
      class ( Current_BSLL_ASC_CSLD_Template ), pointer :: &
        Current_BSLL_ASC_CSLD => null ( )
      type ( StorageDivergenceForm ), allocatable :: &
        StorageDivergence_S
      type ( IncrementDivergence_FV_Form ), dimension ( : ), allocatable :: &
        IncrementDivergence_S
      type ( ApplyDivergence_C_Pointer ) :: &
        ApplyDivergence_S, &
        ApplyDivergence_F
      type ( ApplySources_C_Pointer ) :: &
        ApplySources_S, &
        ApplySources_F
      type ( ApplyRelaxation_C_Pointer ) :: &
        ApplyRelaxation_S, &
        ApplyRelaxation_F
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      DeallocateStorage
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      FinalizeTemplate_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      InitializeTimersStage
    procedure, private, pass :: &
      LoadSolution
    procedure, private, pass :: &
      StoreSolution
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
    procedure, public, pass ( S ) :: &
      LoadSolution_C_BSLL_ASC_CSLD
    procedure, public, pass ( S ) :: &
      StoreSolution_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      InitializeIntermediate_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      IncrementIntermediate_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      ComputeStage_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      IncrementSolution_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      Allocate_RK_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      Deallocate_RK_C_BSLL_ASC_CSLD
  end type Step_RK_C_BSLL_ASC_CSLD_Template

    private :: &
      AllocateStorage

contains


  subroutine InitializeTemplate_C_BSLL_ASC_CSLD &
               ( S, Current_BSLL_ASC_CSLD, NameSuffix, A, B, C )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      Current_BSLL_ASC_CSLD
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iS  !-- iSection

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_BSLL_ASC_CSLD_C_ASC' 

    call S % InitializeTemplate ( NameSuffix, A, B, C )

    S % Bundle_SLL_ASC_CSLD   =>  Current_BSLL_ASC_CSLD % Bundle_SLL_ASC_CSLD
    S % Current_BSLL_ASC_CSLD =>  Current_BSLL_ASC_CSLD

    associate ( B => S % Bundle_SLL_ASC_CSLD )
    S % Chart      =>  B % Base_CSLD
    S % nFibers    =   B % nFibers
    S % nSections  =   B % nSections
    S % Chart_F    =>  B % Fiber_CSLL
    end associate !-- B
    call Show ( S % nFibers,   'nFibers',   S % IGNORABILITY )
    call Show ( S % nSections, 'nSections', S % IGNORABILITY )

    S % Allocated = .false.

    allocate ( S % StorageDivergence_S )
    allocate ( S % IncrementDivergence_S ( S % nSections ) )
    associate ( CB => S % Current_BSLL_ASC_CSLD )
    do iS = 1, S % nSections
      associate ( ID => S % IncrementDivergence_S ( iS ) )
      select type ( CA => CB % Section % Atlas ( iS ) % Element )
      class is ( Current_ASC_Template )
        call ID % Initialize ( CA % Chart )
        call ID % SetStorage ( S % StorageDivergence_S )
      end select !-- CA
      end associate !-- ID
    end do !-- iS
    end associate !-- CB

    allocate ( S % IncrementDamping )
    associate ( ID => S % IncrementDamping )
    call ID % Initialize ( S % Name )
    end associate !-- ID

    S % ApplyDivergence_S % Pointer => ApplyDivergence_C
    S % ApplyDivergence_F % Pointer => null ( )

  end subroutine InitializeTemplate_C_BSLL_ASC_CSLD


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSection

    S % Allocated = .false.

    if ( .not. allocated ( S % IncrementDivergence_S ) ) &
      return

    associate ( ID => S % IncrementDivergence_S ( 1 ) )
    call S % DeallocateMetricDerivatives ( ID )
    end associate !-- ID

    associate ( SD => S % StorageDivergence_S )
    call S % DeallocateStorageDivergence ( SD )
    end associate !-- SD

    do iS = 1, S % nSections
      if ( allocated ( S % BoundaryFluence_CSL_S ) ) then
        associate ( ID => S % IncrementDivergence_S ( iS ) )
        call S % DeallocateBoundaryFluence_CSL &
               ( ID, S % BoundaryFluence_CSL_S ( iS ) % Array )
        end associate !-- ID
      end if
    end do !-- iS

    call S % Deallocate_RK_C_BSLL_ASC_CSLD ( )

  end subroutine DeallocateStorage


  subroutine Compute ( S, Time, TimeStep )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    integer ( KDI ) :: &
      iS     !-- iSection
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_BF

    call Show ( 'Computing a Step_RK_C_BSLL_ASC_CSLD_C_ASC', &
                S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    Timer => PROGRAM_HEADER % Timer ( S % iTimerStep )
    Timer_BF => PROGRAM_HEADER % TimerPointer ( S % iTimerBoundaryFluence )

    if ( associated ( Timer ) ) call Timer % Start ( )

    if ( .not. S % Allocated ) &
      call AllocateStorage ( S )

    if ( associated ( Timer_BF ) ) call Timer_BF % Start ( )
    do iS = 1, S % nSections
      if ( allocated ( S % BoundaryFluence_CSL_S ) ) &
        call S % ClearBoundaryFluence &
               ( S % BoundaryFluence_CSL_S ( iS ) % Array )
    end do !-- iS
    if ( associated ( Timer_BF ) ) call Timer_BF % Stop ( )

    call S % ComputeTemplate ( Time, TimeStep )

    associate ( CB => S % Current_BSLL_ASC_CSLD )

    if ( associated ( Timer_BF ) ) call Timer_BF % Start ( )
    do iS = 1, S % nSections
      select type ( CBA => CB % Section % Atlas ( iS ) % Element )
      class is ( Current_ASC_Template )
        call CBA % AccumulateBoundaryTally &
               ( S % BoundaryFluence_CSL_S ( iS ) % Array )
      end select !-- CBA
    end do !-- iS
    if ( associated ( Timer_BF ) ) call Timer_BF % Stop ( )

    end associate !-- CB, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Compute


  impure elemental subroutine FinalizeTemplate_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    call DeallocateStorage ( S )

    if ( allocated ( S % IncrementDivergence_S ) ) &
      deallocate ( S % IncrementDivergence_S )
    if ( allocated ( S % StorageDivergence_S ) ) &
      deallocate ( S % StorageDivergence_S )

    call S % FinalizeTemplate_C_ASC ( )

  end subroutine FinalizeTemplate_C_BSLL_ASC_CSLD


  subroutine InitializeTimersStage ( S, BaseLevel )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    call PROGRAM_HEADER % AddTimer &
           ( 'Stage', S % iTimerStage, &
             Level = BaseLevel )
      call PROGRAM_HEADER % AddTimer &
             ( 'StoreIntermediate', S % iTimerStoreIntermediate, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'Sections', S % iTimerSections, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'StoreSections', S % iTimerStoreSections, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'Fibers', S % iTimerFibers, Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'LoadSections', S % iTimerLoadSections, &
               Level = BaseLevel + 1 )
      ! call PROGRAM_HEADER % AddTimer &
      !        ( 'ClearIncrement', S % iTimerClearIncrement, &
      !          Level = BaseLevel + 1 )
      ! call PROGRAM_HEADER % AddTimer &
      !        ( 'ApplyDivergence', S % iTimerDivergence, &
      !          Level = BaseLevel + 1 )
      ! call S % InitializeTimersDivergence &
      !        ( BaseLevel + 2 )
      ! call PROGRAM_HEADER % AddTimer &
      !        ( 'ApplySources', S % iTimerSources, &
      !          Level = BaseLevel + 1 )
      ! call PROGRAM_HEADER % AddTimer &
      !        ( 'ApplyRelaxation', S % iTimerRelaxation, &
      !          Level = BaseLevel + 1 )
      ! call PROGRAM_HEADER % AddTimer &
      !        ( 'GhostIncrement', S % iTimerGhost, &
      !          Level = BaseLevel + 1 )

  end subroutine InitializeTimersStage


  subroutine LoadSolution ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    call S % LoadSolution_C_BSLL_ASC_CSLD &
           ( S % Solution_BSLL_ASC_CSLD_S, S % Current_BSLL_ASC_CSLD )

  end subroutine LoadSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    call S % StoreSolution_C_BSLL_ASC_CSLD &
           ( S % Current_BSLL_ASC_CSLD, S % Solution_BSLL_ASC_CSLD_S )

  end subroutine StoreSolution


  subroutine InitializeIntermediate ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    call S % InitializeIntermediate_C_BSLL_ASC_CSLD ( )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    call S % IncrementIntermediate_C_BSLL_ASC_CSLD ( A, iK )

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, Time, TimeStep, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call S % ComputeStage_C_BSLL_ASC_CSLD &
      ( S % Current_BSLL_ASC_CSLD, S % K_BSLL_ASC_CSLD ( iStage ), &
        S % BoundaryFluence_CSL_S, S % Y_BSLL_ASC_CSLD_S, TimeStep, iStage )

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, iS )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    call S % IncrementSolution_C_BSLL_ASC_CSLD ( B, iS )

  end subroutine IncrementSolution


  subroutine LoadSolution_C_BSLL_ASC_CSLD &
               ( Solution_BSLL_ASC_CSLD_S, S, Current_BSLL_ASC_CSLD )

    type ( VariableGroupForm ), dimension ( : ), intent ( inout ) :: &
      Solution_BSLL_ASC_CSLD_S
    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      Current_BSLL_ASC_CSLD

    integer ( KDI ) :: &
      iS     !-- iSection
    class ( CurrentTemplate ), pointer :: &
      Current

    associate ( CB => Current_BSLL_ASC_CSLD )

    do iS = 1, S % nSections
      associate ( Solution => Solution_BSLL_ASC_CSLD_S ( iS ) )
      Current => CB % CurrentSection ( iS )
      call S % LoadSolution_C ( Solution, Current )
      end associate !-- Solution
    end do !-- iS

    end associate !-- CB
    nullify ( Current )

  end subroutine LoadSolution_C_BSLL_ASC_CSLD


  subroutine StoreSolution_C_BSLL_ASC_CSLD &
               ( Current_BSLL_ASC_CSLD, S, Solution_BSLL_ASC_CSLD_S )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      Current_BSLL_ASC_CSLD
    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      Solution_BSLL_ASC_CSLD_S

    integer ( KDI ) :: &
      iS     !-- iSection
    class ( VariableGroupForm ), pointer :: &
      Solution
    class ( CurrentTemplate ), pointer :: &
      Current

    associate ( CB => Current_BSLL_ASC_CSLD )

    do iS = 1, S % nSections
      associate ( Solution => Solution_BSLL_ASC_CSLD_S ( iS ) )
      Current  => CB % CurrentSection ( iS )
      call S % StoreSolution_C ( Current, Solution )
      end associate !-- Solution
    end do !-- iS

    call CB % StoreSections ( )

    end associate !-- CB
    nullify ( Current )

  end subroutine StoreSolution_C_BSLL_ASC_CSLD


  subroutine InitializeIntermediate_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSection

    do iS = 1, S % nSections
      associate &
        ( SV => S % Solution_BSLL_ASC_CSLD_S ( iS ) % Value, &
          YV => S % Y_BSLL_ASC_CSLD_S ( iS ) % Value )
      call Copy ( SV, YV )
      end associate !-- SV, etc.
    end do !-- iS

  end subroutine InitializeIntermediate_C_BSLL_ASC_CSLD


  subroutine IncrementIntermediate_C_BSLL_ASC_CSLD ( S, A, iK )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    integer ( KDI ) :: &
      iS  !-- iSection
    type ( VariableGroupForm ), pointer :: &
      KS

    associate ( KB => S % K_BSLL_ASC_CSLD ( iK ) )

    do iS = 1, S % nSections
      associate ( YS => S % Y_BSLL_ASC_CSLD_S ( iS ) )
      KS => KB % FieldSection ( iS )
      call MultiplyAdd ( YS % Value, KS % Value, A )
      end associate !-- YS
    end do !-- iS

    end associate !-- KB
    nullify ( KS )

  end subroutine IncrementIntermediate_C_BSLL_ASC_CSLD


  subroutine ComputeStage_C_BSLL_ASC_CSLD &
               ( S, C_BSLL_ASC_CSLD, K_BSLL_ASC_CSLD, BF_CSL_S, &
                 Y_BSLL_ASC_CSLD_S, TimeStep, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      C_BSLL_ASC_CSLD
    type ( Storage_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      K_BSLL_ASC_CSLD
    type ( Real_3D_2D_Form ), dimension ( : ), intent ( inout ) :: &
      BF_CSL_S
    type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      Y_BSLL_ASC_CSLD_S
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage
    
    integer ( KDI ) :: &
      iS, &  !-- iSection
      iF     !-- iFiber
    type ( VariableGroupForm ), pointer :: &
      K
    type ( TimerForm ), pointer :: &
      TimerStore, &
      TimerSections, &
      TimerStore_K, &
      TimerFibers, &
      TimerClear, &
      TimerLoad_K
    class ( GeometryFlatForm ), pointer :: &
      G_Base
    class ( CurrentTemplate ), pointer :: &
      C
    type ( IncrementDivergence_FV_Form ) :: &
      ID_Dummy

    associate &
      ( CB => C_BSLL_ASC_CSLD, &
        KB => K_BSLL_ASC_CSLD, &
         B => S % Bundle_SLL_ASC_CSLD )

    if ( iStage > 1 ) then
      TimerStore => PROGRAM_HEADER % TimerPointer &
                      ( S % iTimerStoreIntermediate )
      if ( associated ( TimerStore ) ) call TimerStore % Start ( )
      call S % StoreSolution_C_BSLL_ASC_CSLD &
             ( C_BSLL_ASC_CSLD, Y_BSLL_ASC_CSLD_S )
      if ( associated ( TimerStore ) ) call TimerStore % Stop ( )
    end if

    !-- Sections

    TimerSections => PROGRAM_HEADER % TimerPointer ( S % iTimerSections )
    TimerClear    => PROGRAM_HEADER % TimerPointer ( S % iTimerClearIncrement )

    if ( associated ( TimerSections ) ) call TimerSections % Start ( )

    S % ApplyDivergence_C => S % ApplyDivergence_S % Pointer
    S % ApplySources_C    => S % ApplySources_S    % Pointer
    S % ApplyRelaxation_C => S % ApplyRelaxation_S % Pointer

    do iS = 1, S % nSections

      K => KB % FieldSection ( iS )
      if ( associated ( TimerClear ) ) call TimerClear % Start ( )
      call Clear ( K % Value )
      if ( associated ( TimerClear ) ) call TimerClear % Stop ( )    

      C => CB % CurrentSection ( iS )

      associate ( Chart => S % Chart )

      call S % ComputeStage_C &
             ( S % IncrementDivergence_S ( iS ), C, Chart, K, TimeStep, &
               iStage, GhostExchangeOption = .false. )

      end associate !-- Chart

    end do !-- iS

    S % ApplyRelaxation_C => null ( )
    S % ApplySources_C    => null ( )
    S % ApplyDivergence_C => null ( )

    if ( associated ( TimerSections ) ) call TimerSections % Stop ( )

    !-- Store K to fibers

    TimerStore_K => PROGRAM_HEADER % TimerPointer ( S % iTimerStoreSections )
    if ( associated ( TimerStore_K ) ) call TimerStore_K % Start ( )

    call KB % StoreSections ( )

    if ( associated ( TimerStore_K ) ) call TimerStore_K % Stop ( )

    !-- Fibers

    TimerFibers => PROGRAM_HEADER % TimerPointer ( S % iTimerFibers )
    if ( associated ( TimerFibers ) ) call TimerFibers % Start ( )

    select type ( C_Base => S % Chart )
    class is ( Chart_SL_Template )
      G_Base => C_Base % Geometry ( )
    end select !-- C_Base

    S % ApplyDivergence_C => S % ApplyDivergence_F % Pointer
    S % ApplySources_C    => S % ApplySources_F    % Pointer
    S % ApplyRelaxation_C => S % ApplyRelaxation_F % Pointer

    do iF = 1, S % nFibers

      K => KB % FieldFiber ( iF )

      C => CB % CurrentFiber ( iF )

      associate ( Chart => S % Chart_F )

      call S % ComputeStage_C &
             ( ID_Dummy, C, Chart, K, TimeStep, iStage, &
               GeometryOption = G_Base, &
               iGeometryValueOption = B % iaBaseCell ( iF ) )

      end associate !-- Chart

    end do !-- iF

    S % ApplyRelaxation_C => null ( )
    S % ApplySources_C    => null ( )
    S % ApplyDivergence_C => null ( )

    if ( associated ( TimerFibers ) ) call TimerFibers % Stop ( )

    !-- Load K to sections

    TimerLoad_K => PROGRAM_HEADER % TimerPointer ( S % iTimerLoadSections )
    if ( associated ( TimerLoad_K ) ) call TimerLoad_K % Start ( )

    call KB % LoadSections ( GhostExchangeOption = .true. )

    if ( associated ( TimerLoad_K ) ) call TimerLoad_K % Stop ( )

    end associate !-- CB, etc.
    nullify ( K, C, G_Base )

  end subroutine ComputeStage_C_BSLL_ASC_CSLD


  subroutine IncrementSolution_C_BSLL_ASC_CSLD ( S, B, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iS  !-- iSection
    type ( VariableGroupForm ), pointer :: &
      KS

    associate ( KB => S % K_BSLL_ASC_CSLD ( iStage ) )

    do iS = 1, S % nSections
      associate ( Solution => S % Solution_BSLL_ASC_CSLD_S ( iS ) )
      KS => KB % FieldSection ( iS )
      call MultiplyAdd ( Solution % Value, KS % Value, B )
      end associate !-- Solution
    end do !-- iS

    end associate !-- KB

    nullify ( KS )

  end subroutine IncrementSolution_C_BSLL_ASC_CSLD


  subroutine Allocate_RK_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iStage, &
      iSection
    character ( 1 + 1 ) :: &
      StageNumber
    class ( CurrentTemplate ), pointer :: &
      Current

    allocate ( S % Solution_BSLL_ASC_CSLD_S ( S % nSections ) )
    allocate ( S % Y_BSLL_ASC_CSLD_S ( S % nSections ) )
    allocate ( S % K_BSLL_ASC_CSLD ( S % nStages ) )

    Current => S % Current_BSLL_ASC_CSLD % CurrentSection ( 1 )
    associate &
      ( B => S % Current_BSLL_ASC_CSLD % Bundle_SLL_ASC_CSLD, &
        nEquations => Current % N_CONSERVED, &
        nValues    => Current % nValues )

    do iSection = 1, S % nSections
      call S % Solution_BSLL_ASC_CSLD_S ( iSection ) &
             % Initialize ( [ nValues, nEquations ] )
      call S % Y_BSLL_ASC_CSLD_S ( iSection ) &
             % Initialize ( [ nValues, nEquations ] )
    end do !-- iSection

    do iStage = 1, S % nStages
      write ( StageNumber, fmt = '(a1,i1.1)' ) '_', iStage
      call S % K_BSLL_ASC_CSLD ( iStage ) % Initialize &
             ( B, 'RK_Increment' // StageNumber, nEquations, &
               IgnorabilityOption = CONSOLE % INFO_4 )
    end do !-- iStage

    end associate !-- B, etc.
    nullify ( Current )

  end subroutine Allocate_RK_C_BSLL_ASC_CSLD


  subroutine Deallocate_RK_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % K_BSLL_ASC_CSLD ) ) &
      deallocate ( S % K_BSLL_ASC_CSLD )
    if ( allocated ( S % Y_BSLL_ASC_CSLD_S ) ) &
      deallocate ( S % Y_BSLL_ASC_CSLD_S )
    if ( allocated ( S % Solution_BSLL_ASC_CSLD_S ) ) &
      deallocate ( S % Solution_BSLL_ASC_CSLD_S )
      
  end subroutine Deallocate_RK_C_BSLL_ASC_CSLD


  subroutine AllocateStorage ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSection
    class ( CurrentTemplate ), pointer :: &
      Current_S
    class ( GeometryFlatForm ), pointer :: &
      G

    S % Allocated = .true.

    call S % Allocate_RK_C_BSLL_ASC_CSLD ( )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

      G => Chart % Geometry ( )

      !-- Current_BSLL_ASC_CSLD

      allocate ( S % BoundaryFluence_CSL_S ( S % nSections ) )

      associate ( SD => S % StorageDivergence_S )

      Current_S => S % Current_BSLL_ASC_CSLD % CurrentSection ( 1 )
      call S % AllocateStorageDivergence ( SD, Current_S, G )

      do iS = 1, S % nSections
        associate ( ID => S % IncrementDivergence_S ( iS ) )
        Current_S => S % Current_BSLL_ASC_CSLD % CurrentSection ( iS )
        call S % AllocateBoundaryFluence &
               ( ID, Chart, Current_S % N_CONSERVED, &
                 S % BoundaryFluence_CSL_S ( iS ) % Array )
        end associate !-- ID, etc.
      end do !-- iS

      end associate !-- SD

      nullify ( G )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_BSLL_ASC_CSLD__Template', 'module', &
                  CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    Current_S => S % Current_BSLL_ASC_CSLD % CurrentSection ( 1 )
    call S % AllocateMetricDerivatives &
           ( S % IncrementDivergence_S ( 1 ), Current_S % nValues )

    nullify ( Current_S )

  end subroutine AllocateStorage


end module Step_RK_C_BSLL_ASC_CSLD__Template
