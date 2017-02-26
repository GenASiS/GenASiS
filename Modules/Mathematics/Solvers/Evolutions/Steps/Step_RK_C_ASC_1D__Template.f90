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
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ), abstract :: &
    Step_RK_C_ASC_1D_Template
      integer ( KDI ) :: &
        nCurrents = 0
      type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
        BoundaryFluence_CSL_1D
      logical ( KDL ), dimension ( : ), allocatable :: &
        UseLimiterParameter_1D
      type ( VariableGroupForm ), dimension ( : ), allocatable :: &
        Solution_1D, &
        Y_1D
      type ( VariableGroupForm ), dimension ( :, : ), allocatable :: &
        K_1D
      type ( CurrentPointerForm ), dimension ( : ), allocatable :: &
        Current_1D
      type ( ApplyDivergence_C_Pointer ), dimension ( : ), allocatable :: &
        ApplyDivergence_1D
      type ( ApplySources_C_Pointer ), dimension ( : ), allocatable :: &
        ApplySources_1D
      type ( ApplyRelaxation_C_Pointer ), dimension ( : ), allocatable :: &
        ApplyRelaxation_1D
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_ASC_1D
    procedure, private, pass :: &
      Compute_C_ASC_1D
    generic, public :: &
      Compute => Compute_C_ASC_1D
    procedure, public, pass :: &
      FinalizeTemplate_C_ASC_1D
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
      AllocateStorage, &
      DeallocateStorage

contains


  subroutine InitializeTemplate_C_ASC_1D ( S, NameSuffix, A, B, C, nCurrents )

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
    integer ( KDI ), intent ( in ) :: &
      nCurrents

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_ASC_1D' 

    call S % InitializeTemplate_C_ASC ( NameSuffix, A, B, C )

    S % nCurrents = nCurrents

    allocate ( S % ApplyDivergence_1D ( S % nCurrents ) )
    allocate ( S % ApplySources_1D ( S % nCurrents ) )
    allocate ( S % ApplyRelaxation_1D ( S % nCurrents ) )

  end subroutine InitializeTemplate_C_ASC_1D


  subroutine Compute_C_ASC_1D &
               ( S, Current_ASC_1D, Time, TimeStep, &
                 UseLimiterParameter_1D_Option )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    class ( Current_ASC_ElementForm ), dimension ( : ), intent ( inout ), &
      target :: &
        Current_ASC_1D
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    logical ( KDL ), dimension ( : ), intent ( in ), optional :: &
      UseLimiterParameter_1D_Option

    integer ( KDI ) :: &
      iC  !-- iCurrent

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    allocate ( S % UseLimiterParameter_1D ( S % nCurrents ) )
    S % UseLimiterParameter_1D = .true.
    if ( present ( UseLimiterParameter_1D_Option ) ) &
      S % UseLimiterParameter_1D = UseLimiterParameter_1D_Option

    select type ( Chart => Current_ASC_1D ( 1 ) % Element % Atlas_SC % Chart )
    class is ( Chart_SL_Template )

      S % Grid => Chart

      allocate ( S % Current_1D ( S % nCurrents ) )
      do iC = 1, S % nCurrents
        associate ( CA => Current_ASC_1D ( iC ) % Element )
        S % Current_1D ( iC ) % Pointer => CA % Current ( )
        end associate !-- CA
      end do !-- iC

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC_1D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_C_ASC_1D', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    call AllocateStorage ( S )
    call S % LoadSolution ( S % Solution_1D, S % Current_1D )

    call S % ComputeTemplate ( Time, TimeStep )

    call S % StoreSolution ( S % Current_1D, S % Solution_1D )
    call DeallocateStorage ( S )

    deallocate ( S % Current_1D )
    nullify ( S % Grid )
    deallocate ( S % UseLimiterParameter_1D )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_C_ASC_1D


  impure elemental subroutine FinalizeTemplate_C_ASC_1D ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % ApplyRelaxation_1D ) ) &
      deallocate ( S % ApplyRelaxation_1D )
    if ( allocated ( S % ApplySources_1D ) ) &
      deallocate ( S % ApplySources_1D )
    if ( allocated ( S % ApplyDivergence_1D ) ) &
      deallocate ( S % ApplyDivergence_1D )

    call S % FinalizeTemplate_C_ASC ( )

  end subroutine FinalizeTemplate_C_ASC_1D


  subroutine InitializeIntermediate ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC  !-- iCurrent

    do iC = 1, S % nCurrents
      associate &
        ( SV => S % Solution_1D ( iC ) % Value, &
          YV => S % Y_1D ( iC ) % Value )

      call Copy ( SV, YV )

      end associate !-- SV, etc.
    end do !-- iC

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

    do iC = 1, S % nCurrents
      associate &
        ( YV => S % Y_1D ( iC ) % Value, &
          KV => S % K_1D ( iC, iK ) % Value )

      call MultiplyAdd ( YV, KV, A )

      end associate !-- YV, etc.
    end do !-- iC

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

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeIncrement ) )
    call Timer % Start ( )

    do iC = 1, S % nCurrents
      associate &
        ( C    => S % Current_1D ( iC ) % Pointer, &
          Grid => S % Grid, &
          K    => S % K_1D ( iC, iStage ), &
          BF   => S % BoundaryFluence_CSL_1D ( iC ) % Array, &
          Y    => S % Y_1D ( iC ), &
          ULP  => S % UseLimiterParameter_1D ( iC ) )

      S % ApplyDivergence_C => S % ApplyDivergence_1D ( iC ) % Pointer
      S % ApplySources_C    => S % ApplySources_1D    ( iC ) % Pointer
      S % ApplyRelaxation_C => S % ApplyRelaxation_1D ( iC ) % Pointer

      call S % ComputeStage_C &
             ( C, Grid, K, BF, Y, ULP, TimeStep, iStage )

      S % ApplyRelaxation_C => null ( )
      S % ApplySources_C    => null ( )
      S % ApplyDivergence_C => null ( )

      end associate !-- C, etc.
    end do !-- iC

    call Timer % Stop ( )
    end associate !-- Timer

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

    do iC = 1, S % nCurrents
      associate &
        ( SV => S % Solution_1D ( iC ) % Value, &
          KV => S % K_1D ( iC, iS ) % Value )

      call MultiplyAdd ( SV, KV, B )

      end associate !-- SV, etc.
    end do !-- iC

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

    do iC = 1, S % nStages
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
    character ( LDL ) :: &
      CoordinateSystem

    call S % Allocate_RK_C_1D ( )

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )

      if ( allocated ( S % BoundaryFluence_CSL_1D ) ) &
        deallocate ( S % BoundaryFluence_CSL_1D )
      allocate ( S % BoundaryFluence_CSL_1D ( S % nCurrents ) )

      do iC = 1, S % nCurrents
        call S % AllocateBoundaryFluence &
               ( Grid, S % Current_1D ( iC ) % Pointer % N_CONSERVED, &
                 S % BoundaryFluence_CSL_1D ( iC ) % Array )
      end do !-- iC

      CoordinateSystem = Grid % CoordinateSystem

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    call S % AllocateMetricDerivatives &
           ( CoordinateSystem, S % Current_1D ( 1 ) % Pointer % nValues )

  end subroutine AllocateStorage


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    !-- BoundaryFluence not deallocated here, but instead upon reallocation,
    !   so that its values remain available after Step % Compute

    call S % DeallocateMetricDerivatives ( )
    call S % Deallocate_RK_C_1D ( )

  end subroutine DeallocateStorage


end module Step_RK_C_ASC_1D__Template
