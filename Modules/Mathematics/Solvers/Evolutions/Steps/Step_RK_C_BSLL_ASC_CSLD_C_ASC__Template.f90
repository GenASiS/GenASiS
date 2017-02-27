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
      type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
        BoundaryFluence_CSL_S
      logical ( KDL ) :: &
        UseLimiterParameterSection
      type ( Storage_BSLL_ASC_CSLD_Form ), allocatable :: &
        Solution_BSLL_ASC_CSLD, &
        Y_BSLL_ASC_CSLD
      type ( Storage_BSLL_ASC_CSLD_Form ), dimension ( : ), allocatable :: &
        K_BSLL_ASC_CSLD
      class ( Chart_SLL_Form ), pointer :: &
        Grid_F => null ( )
      class ( Chart_SLD_Form ), pointer :: &
        Grid_S => null ( )
      class ( Current_BSLL_ASC_CSLD_Template ), pointer :: &
        Current_BSLL_ASC_CSLD => null ( )
      type ( ApplyDivergence_C_Pointer ) :: &
        ApplyDivergence_F, &
        ApplyDivergence_S
      type ( ApplySources_C_Pointer ) :: &
        ApplySources_F, &
        ApplySources_S
      type ( ApplyRelaxation_C_Pointer ) :: &
        ApplyRelaxation_F, &
        ApplyRelaxation_S
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC
    procedure, private, pass :: &
      Compute_C_BSLL_ASC_CSLD_C_ASC
    generic, public :: &
      Compute => Compute_C_ASC_1D
    procedure, public, pass :: &
      FinalizeTemplate_C_BSLL_ASC_CSLD_C_ASC
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      IncrementSolution
!    procedure, private, pass :: &
!      ComputeStage
    procedure, private, pass ( S ) :: &
      LoadSolution_C_BSLL_ASC_CSLD
    generic, public :: &
      LoadSolution => LoadSolution_C_BSLL_ASC_CSLD
    procedure, private, pass ( S ) :: &
      StoreSolution_C_BSLL_ASC_CSLD
    generic, public :: &
      StoreSolution => StoreSolution_C_BSLL_ASC_CSLD
    procedure, private, pass ( S ) :: &
      StoreSolution_C_F
    generic, public :: &
      StoreSolution => StoreSolution_C_F
    procedure, public, pass :: &
      InitializeIntermediate_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      IncrementIntermediate_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      IncrementSolution_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      Allocate_RK_C_BSLL_ASC_CSLD
    procedure, public, pass :: &
      Deallocate_RK_C_BSLL_ASC_CSLD
  end type Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template

    private :: &
      AllocateStorage, &
      DeallocateStorage

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

    S % ApplyDivergence_F % Pointer => null ( )

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

    associate ( B => Current_BSLL_ASC_CSLD % Bundle_SLL_ASC_CSLD )
    S % nFibers   =  B % nFibers
    S % nSections =  B % nSections
    S % Grid_F    => B % Fiber_CSLL
    S % Grid_S    => B % Base_CSLD
    end associate !-- B

    S % Current_BSLL_ASC_CSLD => Current_BSLL_ASC_CSLD

    call AllocateStorage ( S )
    call S % LoadSolution ( S % Solution, S % Current )
    call S % LoadSolution ( S % Solution_BSLL_ASC_CSLD, &
                            S % Current_BSLL_ASC_CSLD )

    ! call S % ComputeTemplate ( Time, TimeStep )

    call S % StoreSolution ( S % Current, S % Solution )
    call S % StoreSolution ( S % Current_BSLL_ASC_CSLD, &
                             S % Solution_BSLL_ASC_CSLD )
    call DeallocateStorage ( S )

    S % Current_BSLL_ASC_CSLD => null ( )    
    S % Grid_S => null ( )
    S % Grid_F => null ( )
    S % nSections = 0
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


  subroutine InitializeIntermediate ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    call S % InitializeIntermediate_C_BSLL_ASC_CSLD ( )
    call S % InitializeIntermediate_C ( )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    call S % IncrementIntermediate_C_BSLL_ASC_CSLD ( A, iK )
    call S % IncrementIntermediate_C ( A, iK )

  end subroutine IncrementIntermediate


  subroutine IncrementSolution ( S, B, iS )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    call S % IncrementSolution_C_BSLL_ASC_CSLD ( B, iS )
    call S % IncrementSolution_C ( B, iS )

  end subroutine IncrementSolution


  subroutine LoadSolution_C_BSLL_ASC_CSLD &
               ( Solution_BSLL_ASC_CSLD, S, Current_BSLL_ASC_CSLD )

    type ( Storage_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      Solution_BSLL_ASC_CSLD
    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      Current_BSLL_ASC_CSLD

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iS     !-- iSection
    class ( VariableGroupForm ), pointer :: &
      Solution
    class ( CurrentTemplate ), pointer :: &
      Current

    associate &
      ( CB => Current_BSLL_ASC_CSLD, &
        SB => Solution_BSLL_ASC_CSLD )

    !-- Fibers need Solution when assembling solution from stages
    do iF = 1, S % nFibers
      Current  => CB % CurrentFiber ( iF )
      Solution => SB % FieldFiber ( iF )
      call S % LoadSolution ( Solution, Current )
    end do !-- iC

    !-- Sections need Solution when assembling intermediate state 
    do iS = 1, S % nSections
      Current  => CB % CurrentSection ( iS )
      Solution => SB % FieldSection ( iS )
      call S % LoadSolution ( Solution, Current )
    end do !-- iC

    end associate !-- CB
    nullify ( Current )

  end subroutine LoadSolution_C_BSLL_ASC_CSLD


  subroutine StoreSolution_C_BSLL_ASC_CSLD &
               ( Current_BSLL_ASC_CSLD, S, Solution_BSLL_ASC_CSLD )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      Current_BSLL_ASC_CSLD
    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( in ) :: &
      S
    type ( Storage_BSLL_ASC_CSLD_Form ), intent ( in ) :: &
      Solution_BSLL_ASC_CSLD

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( VariableGroupForm ), pointer :: &
      Solution
    class ( CurrentTemplate ), pointer :: &
      Current

    associate &
      ( CB => Current_BSLL_ASC_CSLD, &
        SB => Solution_BSLL_ASC_CSLD )
    associate ( B => CB % Bundle_SLL_ASC_CSLD )

    !-- New Solution fibers must be stored to Current fibers
    do iF = 1, S % nFibers
      associate ( iBC => B % iaBaseCell ( iF ) )
      Current  => CB % CurrentFiber ( iF )
      Solution => SB % FieldFiber ( iF )
      call S % StoreSolution ( Current, Solution, S % Grid_S, iBC )
      end associate !-- iBC
    end do !-- iC

    !-- Current fibers must be loaded to Current sections
    call CB % LoadSections ( )

    end associate !-- B
    end associate !-- CB
    nullify ( Current )

  end subroutine StoreSolution_C_BSLL_ASC_CSLD


  subroutine StoreSolution_C_F ( Current, S, Solution, Grid_S, iCell_S )

    class ( CurrentTemplate ), intent ( inout ) :: &
      Current
    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( in ) :: &
      Solution
    class ( Chart_SLD_Form ), intent ( in ) :: &
      Grid_S
    integer ( KDI ), intent ( in ) :: &
      iCell_S

    integer ( KDI ) :: &
      iF  !-- iField
    class ( GeometryFlatForm ), pointer :: &
      G_S

    associate ( iaC => Current % iaConserved )
    do iF = 1, Current % N_CONSERVED
      associate &
        ( SV => Solution % Value ( :, iF ), &
          CV => Current % Value ( :, iaC ( iF ) ) )
      call Copy ( SV, CV )
      end associate !-- YV, etc.
    end do !-- iF
    end associate !-- iaC

    G_S => Grid_S % Geometry ( )

    call Current % ComputeFromConserved ( iCell_S, G_S )

    nullify ( G_S )
    
  end subroutine StoreSolution_C_F


  subroutine InitializeIntermediate_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS     !-- iSection
    class ( VariableGroupForm ), pointer :: &
      Solution, &
      Y

    associate &
      ( SB => S % Solution_BSLL_ASC_CSLD, &
        YB => S % Y_BSLL_ASC_CSLD )

    !-- Only need sections for this
    do iS = 1, S % nSections
      Solution => SB % FieldSection ( iS )
      Y => YB % FieldSection ( iS )
      call Copy ( Solution % Value, Y % Value )
    end do !-- iS

    end associate !-- SB, etc.

    nullify ( Solution, Y )

  end subroutine InitializeIntermediate_C_BSLL_ASC_CSLD


  subroutine IncrementIntermediate_C_BSLL_ASC_CSLD ( S, A, iK )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    integer ( KDI ) :: &
      iS     !-- iSection
    class ( VariableGroupForm ), pointer :: &
      Y, &
      K

    associate &
      ( YB => S % Y_BSLL_ASC_CSLD, &
        KB => S % K_BSLL_ASC_CSLD ( iK ) )

    !-- Only need sections for this
    do iS = 1, S % nSections
      Y => YB % FieldSection ( iS )
      K => KB % FieldSection ( iS )
      call MultiplyAdd ( Y % Value, K % Value, A )
    end do !-- iS

    end associate !-- YB, etc.

    nullify ( Y, K )

  end subroutine IncrementIntermediate_C_BSLL_ASC_CSLD


  subroutine IncrementSolution_C_BSLL_ASC_CSLD ( S, B, iS )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iF  !-- iFibers
    class ( VariableGroupForm ), pointer :: &
      Solution, &
      K

    associate &
      ( SB => S % Solution_BSLL_ASC_CSLD, &
        KB => S % K_BSLL_ASC_CSLD ( iS ) )

    !-- Only need fibers for this
    do iF = 1, S % nFibers
      Solution => SB % FieldFiber ( iF )
      K => KB % FieldFiber ( iF )
      call MultiplyAdd ( Solution % Value, K % Value, B )
    end do !-- iF

    end associate !-- SB, etc.

    nullify ( Solution, K )

  end subroutine IncrementSolution_C_BSLL_ASC_CSLD


  subroutine Allocate_RK_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS     !-- iStage
    class ( CurrentTemplate ), pointer :: &
      Current

    allocate ( S % Solution_BSLL_ASC_CSLD )
    allocate ( S % Y_BSLL_ASC_CSLD )
    allocate ( S % K_BSLL_ASC_CSLD ( S % nStages ) )

    Current => S % Current_BSLL_ASC_CSLD % CurrentFiber ( 1 )
    associate &
      ( B => S % Current_BSLL_ASC_CSLD % Bundle_SLL_ASC_CSLD, &
        nEquations => Current % N_CONSERVED )

    call S % Solution_BSLL_ASC_CSLD % Initialize ( B, nEquations )
    call S % Y_BSLL_ASC_CSLD % Initialize ( B, nEquations )
    do iS = 1, S % nStages
      call S % K_BSLL_ASC_CSLD ( iS ) % Initialize ( B, nEquations )
    end do !-- iS

    end associate !-- B, etc.
    nullify ( Current )

  end subroutine Allocate_RK_C_BSLL_ASC_CSLD


  subroutine Deallocate_RK_C_BSLL_ASC_CSLD ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % K_BSLL_ASC_CSLD ) ) &
      deallocate ( S % K_BSLL_ASC_CSLD )
    if ( allocated ( S % Y_BSLL_ASC_CSLD ) ) &
      deallocate ( S % Y_BSLL_ASC_CSLD )
    if ( allocated ( S % Solution_BSLL_ASC_CSLD ) ) &
      deallocate ( S % Solution_BSLL_ASC_CSLD )
      
  end subroutine Deallocate_RK_C_BSLL_ASC_CSLD


  subroutine AllocateStorage ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSection
    character ( LDL ) :: &
      CoordinateSystem
    class ( CurrentTemplate ), pointer :: &
      Current_S

    call S % Allocate_RK_C ( )
    call S % Allocate_RK_C_BSLL_ASC_CSLD ( )

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )

      call S % AllocateBoundaryFluence &
             ( Grid, S % Current % N_CONSERVED, S % BoundaryFluence_CSL )

      if ( allocated ( S % BoundaryFluence_CSL_S ) ) &
        deallocate ( S % BoundaryFluence_CSL_S )
      allocate ( S % BoundaryFluence_CSL_S ( S % nSections ) )

      do iS = 1, S % nSections
        Current_S => S % Current_BSLL_ASC_CSLD % CurrentSection ( iS )
        call S % AllocateBoundaryFluence &
               ( Grid, Current_S % N_CONSERVED, &
                 S % BoundaryFluence_CSL_S ( iS ) % Array )
      end do !-- iC

      CoordinateSystem = Grid % CoordinateSystem

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    call S % AllocateMetricDerivatives &
           ( CoordinateSystem, S % Current % nValues )

    nullify ( Current_S )

  end subroutine AllocateStorage


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S

    !-- BoundaryFluence not deallocated here, but instead upon reallocation,
    !   so that its values remain available after Step % Compute

    call S % DeallocateMetricDerivatives ( )
    call S % Deallocate_RK_C_BSLL_ASC_CSLD ( )
    call S % Deallocate_RK_C ( )

  end subroutine DeallocateStorage


end module Step_RK_C_BSLL_ASC_CSLD_C_ASC__Template
