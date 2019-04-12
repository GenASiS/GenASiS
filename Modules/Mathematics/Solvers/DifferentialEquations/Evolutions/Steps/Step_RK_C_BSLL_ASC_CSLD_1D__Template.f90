!-- Step_RK_C_BSLL_ASC_CSLD_1D is a template for a RungeKutta time step of
!   an array of conserved currents on a bundle.

module Step_RK_C_BSLL_ASC_CSLD_1D__Template

  !-- Step_RungeKutta_Current_BundleSingleLevelDistributed_AtlasSingleChart
  !   _ChartSingleLevelDistributed_1D_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Operations
  use Fields
  use EvolutionBasics
  use Increments
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ), abstract :: &
    Step_RK_C_BSLL_ASC_CSLD_1D_Template
      integer ( KDI ) :: &
        iTimerSections = 0, &
        iTimerStoreSections = 0, &
        iTimerFibers = 0, &
        iTimerLoadSections = 0, &
        nFibers   = 0, &
        nSections = 0, &
        nCurrents = 0
      type ( Real_3D_2D_Form ), dimension ( :, : ), allocatable :: &
        BoundaryFluence_CSL_S
      type ( StorageForm ), dimension ( :, : ), allocatable :: &
        Solution_BSLL_ASC_CSLD_S, &
        Y_BSLL_ASC_CSLD_S
      type ( Storage_BSLL_ASC_CSLD_Form ), dimension ( :, : ), allocatable :: &
        K_BSLL_ASC_CSLD
      class ( Chart_SLL_Form ), pointer :: &
        Chart_F => null ( )
      class ( Bundle_SLL_ASC_CSLD_Form ), pointer :: &
        Bundle_SLL_ASC_CSLD => null ( )
      class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
        pointer :: &
          Current_BSLL_ASC_CSLD_1D => null ( )
      type ( StorageDivergenceForm ), allocatable :: &
        StorageDivergence_S
      type ( IncrementDivergence_FV_Form ), dimension ( :, : ), &
        allocatable :: &
          IncrementDivergence_S
      type ( ComputeConstraints_C_Pointer ), dimension ( : ), allocatable :: &
        ComputeConstraints_S, &
        ComputeConstraints_F
      type ( ApplyDivergence_C_Pointer ), dimension ( : ), allocatable :: &
        ApplyDivergence_S, &
        ApplyDivergence_F
      type ( ApplySources_C_Pointer ), dimension ( : ), allocatable :: &
        ApplySources_S, &
        ApplySources_F
      type ( ApplyRelaxation_C_Pointer ), dimension ( : ), allocatable :: &
        ApplyRelaxation_S, &
        ApplyRelaxation_F
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      DeallocateStorage
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      FinalizeTemplate_C_BSLL_ASC_CSLD_1D
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
      LoadSolution_C_BSLL_ASC_CSLD_1D
    procedure, public, pass ( S ) :: &
      StoreSolution_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      InitializeIntermediate_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      IncrementIntermediate_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      ComputeStage_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      IncrementSolution_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      Allocate_RK_C_BSLL_ASC_CSLD_1D
    procedure, public, pass :: &
      Deallocate_RK_C_BSLL_ASC_CSLD_1D
  end type Step_RK_C_BSLL_ASC_CSLD_1D_Template

    private :: &
      AllocateStorage

contains


  subroutine InitializeTemplate_C_BSLL_ASC_CSLD_1D &
               ( S, I, Current_BSLL_ASC_CSLD_1D, NameSuffix, A, B, C )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    class ( IntegratorHeaderForm ), intent ( in ) :: &
      I
    class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
      intent ( in ), target :: &
        Current_BSLL_ASC_CSLD_1D
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_BSLL_ASC_CSLD_C_ASC_1D' 

    call S % InitializeTemplate ( I, NameSuffix, A, B, C )

    S % Current_BSLL_ASC_CSLD_1D  =>  Current_BSLL_ASC_CSLD_1D
    S % nCurrents = size ( Current_BSLL_ASC_CSLD_1D )
    call Show ( S % nCurrents, 'nCurrents', S % IGNORABILITY )

    S % Bundle_SLL_ASC_CSLD  &
      =>  Current_BSLL_ASC_CSLD_1D ( 1 ) % Element % Bundle_SLL_ASC_CSLD

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
    allocate ( S % IncrementDivergence_S ( S % nSections, S % nCurrents ) )
    do iC = 1, S % nCurrents
      associate ( CB => S % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
      do iS = 1, S % nSections
        associate ( ID => S % IncrementDivergence_S ( iS, iC ) )
        select type ( CA => CB % Section % Atlas ( iS ) % Element )
        class is ( Current_ASC_Template )
          call ID % Initialize ( CA % Chart )
          call ID % SetStorage ( S % StorageDivergence_S )
        end select !-- CA
        end associate !-- ID
      end do !-- iS
      end associate !-- CB
    end do !-- iC

    allocate ( S % ComputeConstraints_S ( S % nCurrents ) )
    allocate ( S % ComputeConstraints_F ( S % nCurrents ) )
    allocate ( S % ApplyDivergence_S ( S % nCurrents ) )
    allocate ( S % ApplyDivergence_F ( S % nCurrents ) )
    allocate ( S % ApplySources_S ( S % nCurrents ) )
    allocate ( S % ApplySources_F ( S % nCurrents ) )
    allocate ( S % ApplyRelaxation_S ( S % nCurrents ) )
    allocate ( S % ApplyRelaxation_F ( S % nCurrents ) )

    do iC = 1, S % nCurrents
      S % ApplyDivergence_S ( iC ) % Pointer => ApplyDivergence_C
      S % ApplyDivergence_F ( iC ) % Pointer => null ( )
    end do !-- iC

  end subroutine InitializeTemplate_C_BSLL_ASC_CSLD_1D


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent

    S % Allocated = .false.

    if ( .not. allocated ( S % IncrementDivergence_S ) ) &
      return

    associate ( ID => S % IncrementDivergence_S ( 1, 1 ) )
    call S % DeallocateMetricDerivatives ( ID )
    end associate !-- ID

    associate ( SD => S % StorageDivergence_S )
    call S % DeallocateStorageDivergence ( SD )
    end associate !-- SD

    do iC = 1, S % nCurrents
      do iS = 1, S % nSections
        if ( allocated ( S % BoundaryFluence_CSL_S ) ) then
          associate ( ID => S % IncrementDivergence_S ( iS, iC ) )
          call S % DeallocateBoundaryFluence_CSL &
                 ( ID, S % BoundaryFluence_CSL_S ( iS, iC ) % Array )
          end associate !-- ID
        end if
      end do !-- iS
    end do !-- iC

    call S % Deallocate_RK_C_BSLL_ASC_CSLD_1D ( )

  end subroutine DeallocateStorage


  subroutine Compute ( S, Time, TimeStep )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_BF

    call Show ( 'Computing a Step_RK_C_BSLL_ASC_CSLD_1D', &
                S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    Timer => PROGRAM_HEADER % Timer ( S % iTimerStep )
    Timer_BF => PROGRAM_HEADER % TimerPointer ( S % iTimerBoundaryFluence )

    if ( associated ( Timer ) ) call Timer % Start ( )

    if ( .not. S % Allocated ) &
      call AllocateStorage ( S )

    if ( associated ( Timer_BF ) ) call Timer_BF % Start ( )
    do iC = 1, S % nCurrents
      do iS = 1, S % nSections
        if ( allocated ( S % BoundaryFluence_CSL_S ) ) &
          call S % ClearBoundaryFluence &
                 ( S % BoundaryFluence_CSL_S ( iS, iC ) % Array )
      end do !-- iS
    end do !-- iC
    if ( associated ( Timer_BF ) ) call Timer_BF % Stop ( )

    call S % ComputeTemplate ( Time, TimeStep )

    if ( associated ( Timer_BF ) ) call Timer_BF % Start ( )
    do iC = 1, S % nCurrents
      associate ( CB => S % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
      do iS = 1, S % nSections
        select type ( CBA => CB % Section % Atlas ( iS ) % Element )
        class is ( Current_ASC_Template )
          call CBA % AccumulateBoundaryTally &
                 ( S % BoundaryFluence_CSL_S ( iS, iC ) % Array )
        end select !-- CBA
      end do !-- iS
      end associate !-- CB, etc.
    end do !-- iC
    if ( associated ( Timer_BF ) ) call Timer_BF % Stop ( )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Compute


  impure elemental subroutine FinalizeTemplate_C_BSLL_ASC_CSLD_1D ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    call DeallocateStorage ( S )

    if ( allocated ( S % ApplyRelaxation_F ) ) &
      deallocate ( S % ApplyRelaxation_F )
    if ( allocated ( S % ApplyRelaxation_S ) ) &
      deallocate ( S % ApplyRelaxation_S )
    if ( allocated ( S % ApplySources_F ) ) &
      deallocate ( S % ApplySources_F )
    if ( allocated ( S % ApplySources_S ) ) &
      deallocate ( S % ApplySources_S )
    if ( allocated ( S % ApplyDivergence_F ) ) &
      deallocate ( S % ApplyDivergence_F )
    if ( allocated ( S % ApplyDivergence_S ) ) &
      deallocate ( S % ApplyDivergence_S )
    if ( allocated ( S % ComputeConstraints_F ) ) &
      deallocate ( S % ComputeConstraints_F )
    if ( allocated ( S % ComputeConstraints_S ) ) &
      deallocate ( S % ComputeConstraints_S )
    if ( allocated ( S % IncrementDivergence_S ) ) &
      deallocate ( S % IncrementDivergence_S )
    if ( allocated ( S % StorageDivergence_S ) ) &
      deallocate ( S % StorageDivergence_S )

    nullify ( S % Current_BSLL_ASC_CSLD_1D )

    call S % FinalizeTemplate_C_ASC ( )

  end subroutine FinalizeTemplate_C_BSLL_ASC_CSLD_1D


  subroutine InitializeTimersStage ( S, BaseLevel )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
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

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    call S % LoadSolution_C_BSLL_ASC_CSLD_1D &
           ( S % Solution_BSLL_ASC_CSLD_S, S % Current_BSLL_ASC_CSLD_1D )

  end subroutine LoadSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    call S % StoreSolution_C_BSLL_ASC_CSLD_1D &
           ( S % Current_BSLL_ASC_CSLD_1D, S % Solution_BSLL_ASC_CSLD_S )

  end subroutine StoreSolution


  subroutine InitializeIntermediate ( S, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iStage

    if ( iStage > 1 ) &
      call S % InitializeIntermediate_C_BSLL_ASC_CSLD_1D ( )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    call S % IncrementIntermediate_C_BSLL_ASC_CSLD_1D ( A, iK )

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, Time, TimeStep, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call S % ComputeStage_C_BSLL_ASC_CSLD_1D &
      ( S % Current_BSLL_ASC_CSLD_1D, S % K_BSLL_ASC_CSLD ( :, iStage ), &
        S % BoundaryFluence_CSL_S, S % Y_BSLL_ASC_CSLD_S, TimeStep, iStage )

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, iS )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    call S % IncrementSolution_C_BSLL_ASC_CSLD_1D ( B, iS )

  end subroutine IncrementSolution


  subroutine LoadSolution_C_BSLL_ASC_CSLD_1D &
               ( Solution_BSLL_ASC_CSLD_S, S, Current_BSLL_ASC_CSLD_1D )

    type ( StorageForm ), dimension ( :, : ), intent ( inout ) :: &
      Solution_BSLL_ASC_CSLD_S
    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( in ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
      intent ( in ) :: &
        Current_BSLL_ASC_CSLD_1D

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent
    class ( CurrentTemplate ), pointer :: &
      Current

    do iC = 1, S % nCurrents
      associate ( CB => Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
      do iS = 1, S % nSections
        associate ( Solution => Solution_BSLL_ASC_CSLD_S ( iS, iC ) )
        Current => CB % CurrentSection ( iS )
        call S % LoadSolution_C ( Solution, Current )
        end associate !-- Solution
      end do !-- iS
      end associate !-- CB
    end do !-- iC
    nullify ( Current )

  end subroutine LoadSolution_C_BSLL_ASC_CSLD_1D


  subroutine StoreSolution_C_BSLL_ASC_CSLD_1D &
               ( Current_BSLL_ASC_CSLD_1D, S, Solution_BSLL_ASC_CSLD_S )

    class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
      intent ( inout ) :: &
        Current_BSLL_ASC_CSLD_1D
    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( in ) :: &
      S
    type ( StorageForm ), dimension ( :, : ), intent ( in ) :: &
      Solution_BSLL_ASC_CSLD_S

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent
    class ( CurrentTemplate ), pointer :: &
      Current

    do iC = 1, S % nCurrents
      associate ( CB => Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
      do iS = 1, S % nSections
        associate ( Solution => Solution_BSLL_ASC_CSLD_S ( iS, iC ) )
        Current  => CB % CurrentSection ( iS )
        call S % StoreSolution_C ( Current, Solution )
        end associate !-- Solution
      end do !-- iS
      call CB % StoreSections ( )
      end associate !-- CB
    end do !-- iC
    nullify ( Current )

  end subroutine StoreSolution_C_BSLL_ASC_CSLD_1D


  subroutine InitializeIntermediate_C_BSLL_ASC_CSLD_1D ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent

    do iC = 1, S % nCurrents
      do iS = 1, S % nSections
        associate &
          ( SV => S % Solution_BSLL_ASC_CSLD_S ( iS, iC ) % Value, &
            YV => S % Y_BSLL_ASC_CSLD_S ( iS, iC ) % Value )
        call Copy ( SV, YV )
        end associate !-- SV, etc.
      end do !-- iS
    end do !-- iC

  end subroutine InitializeIntermediate_C_BSLL_ASC_CSLD_1D


  subroutine IncrementIntermediate_C_BSLL_ASC_CSLD_1D ( S, A, iK )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent
    type ( StorageForm ), pointer :: &
      KS

    do iC = 1, S % nCurrents
      associate ( KB => S % K_BSLL_ASC_CSLD ( iC, iK ) )
      do iS = 1, S % nSections
        associate ( YS => S % Y_BSLL_ASC_CSLD_S ( iS, iC ) )
        KS => KB % FieldSection ( iS )
        call MultiplyAdd ( YS % Value, KS % Value, A )
        end associate !-- YS
      end do !-- iS
      end associate !-- KB
    end do !-- iC
    nullify ( KS )

  end subroutine IncrementIntermediate_C_BSLL_ASC_CSLD_1D


  subroutine ComputeStage_C_BSLL_ASC_CSLD_1D &
               ( S, C_BSLL_ASC_CSLD_1D, K_BSLL_ASC_CSLD, BF_CSL_S, &
                 Y_BSLL_ASC_CSLD_S, TimeStep, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
      intent ( inout ) :: &
        C_BSLL_ASC_CSLD_1D
    type ( Storage_BSLL_ASC_CSLD_Form ), dimension ( : ), intent ( inout ) :: &
      K_BSLL_ASC_CSLD
    type ( Real_3D_2D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF_CSL_S
    type ( StorageForm ), dimension ( :, : ), intent ( in ) :: &
      Y_BSLL_ASC_CSLD_S
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage
    
    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iS, &  !-- iSection
      iC     !-- iCurrent
    type ( StorageForm ), pointer :: &
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

    if ( iStage > 1 ) then
      TimerStore => PROGRAM_HEADER % TimerPointer &
                      ( S % iTimerStoreIntermediate )
      if ( associated ( TimerStore ) ) call TimerStore % Start ( )
      call S % StoreSolution_C_BSLL_ASC_CSLD_1D &
             ( C_BSLL_ASC_CSLD_1D, Y_BSLL_ASC_CSLD_S )
      if ( associated ( TimerStore ) ) call TimerStore % Stop ( )
    end if

    !-- Sections

    TimerSections &
      => PROGRAM_HEADER % TimerPointer ( S % iTimerSections )
    TimerClear &
      => PROGRAM_HEADER % TimerPointer ( S % iTimerClearIncrement )

    if ( associated ( TimerSections ) ) call TimerSections % Start ( )

    do iC = 1, S % nCurrents
      associate &
        ( CB => C_BSLL_ASC_CSLD_1D ( iC ) % Element, &
          KB => K_BSLL_ASC_CSLD ( iC ) )

      S % ApplyDivergence_C => S % ApplyDivergence_S ( iC ) % Pointer
      S % ApplySources_C    => S % ApplySources_S ( iC )    % Pointer
      S % ApplyRelaxation_C => S % ApplyRelaxation_S ( iC ) % Pointer

      do iS = 1, S % nSections
 
        K => KB % FieldSection ( iS )
        if ( associated ( TimerClear ) ) call TimerClear % Start ( )
        call Clear ( K % Value )
        if ( associated ( TimerClear ) ) call TimerClear % Stop ( )    

        C => CB % CurrentSection ( iS )

        associate ( Chart => S % Chart )

        call S % ComputeStage_C &
               ( S % IncrementDivergence_S ( iS, iC ), C, Chart, K, TimeStep, &
                 iStage, GhostExchangeOption = .false. )

        end associate !-- Chart

      end do !-- iS

      S % ApplyRelaxation_C => null ( )
      S % ApplySources_C    => null ( )
      S % ApplyDivergence_C => null ( )

      if ( associated ( TimerSections ) ) call TimerSections % Stop ( )

      end associate !-- CB, etc.

    end do !-- iC

    !-- Store K to fibers

    TimerStore_K => PROGRAM_HEADER % TimerPointer ( S % iTimerStoreSections )
    if ( associated ( TimerStore_K ) ) call TimerStore_K % Start ( )

    do iC = 1, S % nCurrents
      associate ( KB => K_BSLL_ASC_CSLD ( iC ) )
      call KB % StoreSections ( )
      end associate !-- KB
    end do !-- iC

    if ( associated ( TimerStore_K ) ) call TimerStore_K % Stop ( )

    !-- Fibers

    TimerFibers => PROGRAM_HEADER % TimerPointer ( S % iTimerFibers )
    if ( associated ( TimerFibers ) ) call TimerFibers % Start ( )

    select type ( C_Base => S % Chart )
    class is ( Chart_SL_Template )
      G_Base => C_Base % Geometry ( )
    end select !-- C_Base

    do iC = 1, S % nCurrents
      associate &
        ( CB => C_BSLL_ASC_CSLD_1D ( iC ) % Element, &
          KB => K_BSLL_ASC_CSLD ( iC ), &
           B => S % Bundle_SLL_ASC_CSLD )

      S % ApplyDivergence_C => S % ApplyDivergence_F ( iC ) % Pointer
      S % ApplySources_C    => S % ApplySources_F ( iC )    % Pointer
      S % ApplyRelaxation_C => S % ApplyRelaxation_F ( iC ) % Pointer

      do iF = 1, S % nFibers

        K => KB % FieldFiber ( iF )

        C => CB % CurrentFiber ( iF )

        associate ( Chart => S % Chart_F )

        call S % ComputeStage_C &
               ( ID_Dummy, C, Chart, K, TimeStep, iStage, &
                 GeometryOption = G_Base, &
                 iStrgeometryValueOption = B % iaBaseCell ( iF ) )

        end associate !-- Chart

      end do !-- iF

      S % ApplyRelaxation_C => null ( )
      S % ApplySources_C    => null ( )
      S % ApplyDivergence_C => null ( )

      end associate !-- CB, etc.
    end do !-- iC

    if ( associated ( TimerFibers ) ) call TimerFibers % Stop ( )

    !-- Load K to sections

    TimerLoad_K => PROGRAM_HEADER % TimerPointer ( S % iTimerLoadSections )
    if ( associated ( TimerLoad_K ) ) call TimerLoad_K % Start ( )

    do iC = 1, S % nCurrents
      associate ( KB => K_BSLL_ASC_CSLD ( iC ) )
      call KB % LoadSections ( GhostExchangeOption = .true. )
      end associate !-- KB
    end do !-- iC

    if ( associated ( TimerLoad_K ) ) call TimerLoad_K % Stop ( )

    nullify ( K, C, G_Base )

  end subroutine ComputeStage_C_BSLL_ASC_CSLD_1D


  subroutine IncrementSolution_C_BSLL_ASC_CSLD_1D ( S, B, iStage )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent
    type ( StorageForm ), pointer :: &
      KS

    do iC = 1, S % nCurrents
      associate ( KB => S % K_BSLL_ASC_CSLD ( iC, iStage ) )
      do iS = 1, S % nSections
        associate ( Solution => S % Solution_BSLL_ASC_CSLD_S ( iS, iC ) )
        KS => KB % FieldSection ( iS )
        call MultiplyAdd ( Solution % Value, KS % Value, B )
        end associate !-- Solution
      end do !-- iS
      end associate !-- KB
    end do !-- iC

    nullify ( KS )

  end subroutine IncrementSolution_C_BSLL_ASC_CSLD_1D


  subroutine Allocate_RK_C_BSLL_ASC_CSLD_1D ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iStage, &
      iSection, &
      iC  !-- iCurrent
    character ( 1 + 1 ) :: &
      StageNumber
    class ( CurrentTemplate ), pointer :: &
      Current

    allocate ( S % Solution_BSLL_ASC_CSLD_S ( S % nSections, S % nCurrents ) )
    allocate ( S % Y_BSLL_ASC_CSLD_S ( S % nSections, S % nCurrents ) )
    allocate ( S % K_BSLL_ASC_CSLD ( S % nCurrents, S % nStages ) )

    Current => S % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element &
                 % CurrentSection ( 1 )
    associate &
      ( B => S % Bundle_SLL_ASC_CSLD, &
        nEquations => Current % N_CONSERVED, &
        nValues    => Current % nValues )

    do iC = 1, S % nCurrents

      do iSection = 1, S % nSections
        call S % Solution_BSLL_ASC_CSLD_S ( iSection, iC ) &
               % Initialize ( [ nValues, nEquations ] )
        call S % Y_BSLL_ASC_CSLD_S ( iSection, iC ) &
               % Initialize ( [ nValues, nEquations ] )
      end do !-- iSection

      do iStage = 1, S % nStages
        write ( StageNumber, fmt = '(a1,i1.1)' ) '_', iStage
        call S % K_BSLL_ASC_CSLD ( iC, iStage ) % Initialize &
               ( B, 'RK_Increment' // StageNumber, nEquations, &
                 IgnorabilityOption = CONSOLE % INFO_4 )
      end do !-- iStage

    end do !-- iC

    end associate !-- B, etc.
    nullify ( Current )

  end subroutine Allocate_RK_C_BSLL_ASC_CSLD_1D


  subroutine Deallocate_RK_C_BSLL_ASC_CSLD_1D ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % K_BSLL_ASC_CSLD ) ) &
      deallocate ( S % K_BSLL_ASC_CSLD )
    if ( allocated ( S % Y_BSLL_ASC_CSLD_S ) ) &
      deallocate ( S % Y_BSLL_ASC_CSLD_S )
    if ( allocated ( S % Solution_BSLL_ASC_CSLD_S ) ) &
      deallocate ( S % Solution_BSLL_ASC_CSLD_S )
      
  end subroutine Deallocate_RK_C_BSLL_ASC_CSLD_1D


  subroutine AllocateStorage ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent
    class ( CurrentTemplate ), pointer :: &
      Current_S
    class ( GeometryFlatForm ), pointer :: &
      G

    S % Allocated = .true.

    call S % Allocate_RK_C_BSLL_ASC_CSLD_1D ( )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

      G => Chart % Geometry ( )

      !-- Current_BSLL_ASC_CSLD

      associate ( SD => S % StorageDivergence_S )
      Current_S => S % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element &
                     % CurrentSection ( 1 )
      call S % AllocateStorageDivergence ( SD, Current_S, G )
      end associate !-- SD

      allocate ( S % BoundaryFluence_CSL_S ( S % nSections, S % nCurrents ) )
      do iC = 1, S % nCurrents
        do iS = 1, S % nSections
          associate ( ID => S % IncrementDivergence_S ( iS, iC ) )
          Current_S => S % Current_BSLL_ASC_CSLD_1D ( iC ) % Element &
                         % CurrentSection ( iS )
          call S % AllocateBoundaryFluence &
                 ( ID, Current_S, Chart, &
                   S % BoundaryFluence_CSL_S ( iS, iC ) % Array )
          end associate !-- ID, etc.
        end do !-- iS
      end do !-- iC

      nullify ( G )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_BSLL_ASC_CSLD_1D__Template', 'module', &
                  CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    Current_S => S % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element &
                   % CurrentSection ( 1 )
    call S % AllocateMetricDerivatives &
           ( S % IncrementDivergence_S ( 1, 1 ), Current_S )

    nullify ( Current_S )

  end subroutine AllocateStorage


end module Step_RK_C_BSLL_ASC_CSLD_1D__Template
