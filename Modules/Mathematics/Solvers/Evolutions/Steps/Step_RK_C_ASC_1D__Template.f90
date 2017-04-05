!-- Step_RK_C_ASC is a template for a RungeKutta time step of an array of
!   conserved currents.

module Step_RK_C_ASC_1D__Template

  !-- Step_RungeKutta_Current_AtlasSingleChart_1D_Template

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
    Step_RK_C_ASC_1D_Template
      integer ( KDI ) :: &
        nCurrents = 0
      type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
        BoundaryFluence_CSL_1D
      type ( VariableGroupForm ), dimension ( : ), allocatable :: &
        Solution_1D, &
        Y_1D
      type ( VariableGroupForm ), dimension ( :, : ), allocatable :: &
        K_1D
      type ( CurrentPointerForm ), dimension ( : ), allocatable :: &
        Current_1D
      class ( Current_ASC_ElementForm ), dimension ( : ), pointer :: &
        Current_ASC_1D => null ( )
      type ( StorageDivergenceForm ), dimension ( : ), allocatable :: &
        StorageDivergence_1D
      type ( IncrementDivergence_FV_Form ), dimension ( : ), allocatable :: &
        IncrementDivergence_1D
      type ( ApplyDivergence_C_Pointer ), dimension ( : ), allocatable :: &
        ApplyDivergence_1D
      type ( ApplySources_C_Pointer ), dimension ( : ), allocatable :: &
        ApplySources_1D
      type ( ApplyRelaxation_C_Pointer ), dimension ( : ), allocatable :: &
        ApplyRelaxation_1D
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_ASC_1D
    procedure, public, pass :: &
      DeallocateStorage
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      FinalizeTemplate_C_ASC_1D
    procedure, private, pass :: &
      InitializeTimersDivergence
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
    procedure, private, pass ( S ) :: &
      LoadSolution_C_1D
    generic, public :: &
      LoadSolution => LoadSolution_C_1D
    procedure, private, pass ( S ) :: &
      StoreSolution_C_1D
    generic, public :: &
      StoreSolution => StoreSolution_C_1D
    procedure, public, pass :: &
      Allocate_RK_C_1D
    procedure, public, pass :: &
      Deallocate_RK_C_1D
  end type Step_RK_C_ASC_1D_Template

    private :: &
      AllocateStorage

contains


  subroutine InitializeTemplate_C_ASC_1D &
               ( S, Current_ASC_1D, NameSuffix, A, B, C )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    class ( Current_ASC_ElementForm ), dimension ( : ), intent ( in ), &
      target :: &
        Current_ASC_1D
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iC  !-- iCurrent

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_ASC_1D' 

    call S % InitializeTemplate ( NameSuffix, A, B, C )

    S % Current_ASC_1D => Current_ASC_1D
    S % nCurrents = size ( Current_ASC_1D )
    call Show ( S % nCurrents, 'nCurrents', S % IGNORABILITY )

    select type ( Chart => Current_ASC_1D ( 1 ) % Element % Atlas_SC % Chart )
    class is ( Chart_SL_Template )

      S % Chart => Chart

      allocate ( S % Current_1D ( S % nCurrents ) )
      do iC = 1, S % nCurrents
        associate ( CA => Current_ASC_1D ( iC ) % Element )
        S % Current_1D ( iC ) % Pointer => CA % Current ( )
        end associate !-- CA
      end do !-- iC

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC_1D__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C_ASC_1D', 'subroutine', &
                  CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    S % Allocated = .false.

    allocate ( S % StorageDivergence_1D ( S % nCurrents ) )
    allocate ( S % IncrementDivergence_1D ( S % nCurrents ) )
    do iC = 1, S % nCurrents
      associate ( ID => S % IncrementDivergence_1D ( iC ) )
      call ID % Initialize ( Current_ASC_1D ( iC ) % Element % Chart )
      call ID % SetStorage ( S % StorageDivergence_1D ( iC ) )
      end associate !-- ID
    end do !-- iC

    allocate ( S % IncrementDamping )
    associate ( ID => S % IncrementDamping )
    call ID % Initialize ( S % Name )
    end associate !-- ID

    allocate ( S % ApplyDivergence_1D ( S % nCurrents ) )
    allocate ( S % ApplySources_1D ( S % nCurrents ) )
    allocate ( S % ApplyRelaxation_1D ( S % nCurrents ) )

  end subroutine InitializeTemplate_C_ASC_1D


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC  !-- iCurrent

    S % Allocated = .false.

    if ( .not. allocated ( S % IncrementDivergence_1D ) ) &
      return

    call S % DeallocateMetricDerivatives ( S % IncrementDivergence_1D ( 1 ) )

    if ( .not. allocated ( S % BoundaryFluence_CSL_1D ) ) &
      return

    do iC = 1, S % nCurrents
      associate &
        ( ID => S % IncrementDivergence_1D ( iC ), &
          SD => S % StorageDivergence_1D ( iC ) )
      call S % DeallocateBoundaryFluence_CSL &
             ( ID, S % BoundaryFluence_CSL_1D ( iC ) % Array )
      call S % DeallocateStorageDivergence ( SD )
      end associate !-- ID
    end do !-- iC

    call S % Deallocate_RK_C_1D ( )

  end subroutine DeallocateStorage


  subroutine Compute ( S, Time, TimeStep )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % Timer ( S % iTimerStep )
    if ( associated ( Timer ) ) call Timer % Start ( )

    if ( .not. S % Allocated ) &
      call AllocateStorage ( S )

    do iC = 1, S % nCurrents
      if ( allocated ( S % BoundaryFluence_CSL_1D ) ) &
        call S % ClearBoundaryFluence &
               ( S % BoundaryFluence_CSL_1D ( iC ) % Array )
    end do !-- iC

    call S % LoadSolution ( S % Solution_1D, S % Current_1D )
    call S % ComputeTemplate ( Time, TimeStep )
    call S % StoreSolution ( S % Current_1D, S % Solution_1D )

    do iC = 1, S % nCurrents
      associate ( CA => S % Current_ASC_1D ( iC ) % Element )
      if ( allocated ( S % BoundaryFluence_CSL_1D ) ) &
        call CA % AccumulateBoundaryTally &
               ( S % BoundaryFluence_CSL_1D ( iC ) % Array )
      end associate !-- CA
    end do !-- iC

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Compute


  impure elemental subroutine FinalizeTemplate_C_ASC_1D ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    call DeallocateStorage ( S )

    if ( allocated ( S % ApplyRelaxation_1D ) ) &
      deallocate ( S % ApplyRelaxation_1D )
    if ( allocated ( S % ApplySources_1D ) ) &
      deallocate ( S % ApplySources_1D )
    if ( allocated ( S % ApplyDivergence_1D ) ) &
      deallocate ( S % ApplyDivergence_1D )
    if ( allocated ( S % IncrementDivergence_1D ) ) &
      deallocate ( S % IncrementDivergence_1D )
    if ( allocated ( S % StorageDivergence_1D ) ) &
      deallocate ( S % StorageDivergence_1D )
    if ( allocated ( S % Current_1D ) ) &
      deallocate ( S % Current_1D )

    nullify ( S % Current_ASC_1D )

    call S % FinalizeTemplate_C_ASC ( )

  end subroutine FinalizeTemplate_C_ASC_1D


  subroutine InitializeTimersDivergence ( S, BaseLevel )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    integer ( KDI ) :: &
      iC  !-- iCurrent

    do iC = 1, S % nCurrents
      associate &
        ( SD => S % StorageDivergence_1D ( iC ), &
          C  => S % Current_1D ( iC ) % Pointer )
      call SD % InitializeTimers ( C % Name, BaseLevel )
      end associate !-- SD, etc.
    end do

  end subroutine InitializeTimersDivergence


  subroutine InitializeIntermediate ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer &
               ( S % iTimerInitializeIntermediate )
    if ( associated ( Timer ) ) call Timer % Start ( )

    do iC = 1, S % nCurrents
      associate &
        ( SV => S % Solution_1D ( iC ) % Value, &
          YV => S % Y_1D ( iC ) % Value )
      call Copy ( SV, YV )
      end associate !-- SV, etc.
    end do !-- iC

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerIncrementIntermediate )
    if ( associated ( Timer ) ) call Timer % Start ( )

    do iC = 1, S % nCurrents
      associate &
        ( YV => S % Y_1D ( iC ) % Value, &
          KV => S % K_1D ( iC, iK ) % Value )
      call MultiplyAdd ( YV, KV, A )
      end associate !-- YV, etc.
    end do !-- iC

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, Time, TimeStep, iStage )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerStage )
    if ( associated ( Timer ) ) call Timer % Start ( )

    do iC = 1, S % nCurrents
      associate &
        ( C     => S % Current_1D ( iC ) % Pointer, &
          Chart => S % Chart, &
          K     => S % K_1D ( iC, iStage ), &
          Y     => S % Y_1D ( iC ) )

      if ( iStage > 1 ) &
        call S % StoreSolution ( C, Y )

      call Clear ( K % Value )

      S % ApplyDivergence_C => S % ApplyDivergence_1D ( iC ) % Pointer
      S % ApplySources_C    => S % ApplySources_1D    ( iC ) % Pointer
      S % ApplyRelaxation_C => S % ApplyRelaxation_1D ( iC ) % Pointer

      call S % ComputeStage_C &
             ( S % IncrementDivergence_1D ( iC ), C, Chart, K, TimeStep, &
               iStage )

      S % ApplyRelaxation_C => null ( )
      S % ApplySources_C    => null ( )
      S % ApplyDivergence_C => null ( )

      end associate !-- C, etc.
    end do !-- iC

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, iS )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerIncrementSolution )
    if ( associated ( Timer ) ) call Timer % Start ( )

    do iC = 1, S % nCurrents
      associate &
        ( SV => S % Solution_1D ( iC ) % Value, &
          KV => S % K_1D ( iC, iS ) % Value )
      call MultiplyAdd ( SV, KV, B )
      end associate !-- SV, etc.
    end do !-- iC

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine IncrementSolution


  subroutine LoadSolution_C_1D ( Solution_1D, S, Current_1D )

    type ( VariableGroupForm ), dimension ( : ), intent ( inout ) :: &
      Solution_1D
    class ( Step_RK_C_ASC_1D_Template ), intent ( in ) :: &
      S
    type ( CurrentPointerForm ), dimension ( : ), intent ( in ) :: &
      Current_1D

    integer ( KDI ) :: &
      iC  !-- iCurrent

    do iC = 1, size ( Current_1D )
      call S % LoadSolution &
             ( Solution_1D ( iC ), Current_1D ( iC ) % Pointer )
    end do !-- iC

  end subroutine LoadSolution_C_1D


  subroutine StoreSolution_C_1D ( Current_1D, S, Solution_1D )

    type ( CurrentPointerForm ), dimension ( : ), intent ( in ) :: &
      Current_1D
    class ( Step_RK_C_ASC_1D_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), dimension ( : ), intent ( inout ) :: &
      Solution_1D

    integer ( KDI ) :: &
      iC  !-- iCurrent

    do iC = 1, size ( Current_1D )
      call S % StoreSolution &
             ( Current_1D ( iC ) % Pointer, Solution_1D ( iC ) )
    end do !-- iC

  end subroutine StoreSolution_C_1D


  subroutine Allocate_RK_C_1D ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC, &  !-- iCurrent
      iS     !-- iStage


    allocate ( S % Solution_1D ( S % nCurrents ) )
    allocate ( S % Y_1D ( S % nCurrents ) )
    allocate ( S % K_1D ( S % nCurrents, S % nStages ) )

    do iC = 1, S % nCurrents
      associate &
        ( nEquations => S % Current_1D ( iC ) % Pointer % N_CONSERVED, &
          nValues    => S % Current_1D ( iC ) % Pointer % nValues )

      call S % Solution_1D ( iC ) % Initialize ( [ nValues, nEquations ] )
      call S % Y_1D ( iC ) % Initialize ( [ nValues, nEquations ] )
      do iS = 1, S % nStages
        call S % K_1D ( iC, iS ) % Initialize ( [ nValues, nEquations ] )
      end do !-- iS

      end associate !-- nEquations, etc.
    end do

  end subroutine Allocate_RK_C_1D


  subroutine Deallocate_RK_C_1D ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % K_1D ) ) &
      deallocate ( S % K_1D )
    if ( allocated ( S % Y_1D ) ) &
      deallocate ( S % Y_1D )
    if ( allocated ( S % Solution_1D ) ) &
      deallocate ( S % Solution_1D )
      
  end subroutine Deallocate_RK_C_1D


  subroutine AllocateStorage ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC
    class ( GeometryFlatForm ), pointer :: &
      G

    S % Allocated = .true.

    call S % Allocate_RK_C_1D ( )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

      G => Chart % Geometry ( )

      allocate ( S % BoundaryFluence_CSL_1D ( S % nCurrents ) )

      do iC = 1, S % nCurrents
        associate &
          ( ID => S % IncrementDivergence_1D ( iC ), &
            SD => S % StorageDivergence_1D ( iC ), &
            C  => S % Current_1D ( iC ) % Pointer )

        call S % AllocateStorageDivergence ( SD, C, G )

        call S % AllocateBoundaryFluence &
               ( ID, Chart, C % N_CONSERVED, &
                 S % BoundaryFluence_CSL_1D ( iC ) % Array )

        end associate !-- ID, etc.
      end do !-- iC

      call S % AllocateMetricDerivatives &
             ( S % IncrementDivergence_1D ( 1 ), &
               S % Current_1D ( 1 ) % Pointer % nValues )

      nullify ( G )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

  end subroutine AllocateStorage


end module Step_RK_C_ASC_1D__Template
