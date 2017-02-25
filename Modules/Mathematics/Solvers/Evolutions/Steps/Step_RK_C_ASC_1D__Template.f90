!-- Step_RK_C_ASC is a template for a RungeKutta time step of an array of
!   conserved currents.

module Step_RK_C_ASC_1D__Template

  !-- Step_RungeKutta_Current_AtlasSingleChart_1D_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Fields
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ), abstract :: &
    Step_RK_C_ASC_1D_Template
      integer ( KDI ) :: &
        nCurrents = 0
      type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
        BoundaryFluence_CSL_C_1D
      logical ( KDL ), dimension ( : ), allocatable :: &
        UseLimiterParameter_C_1D
      type ( VariableGroupForm ), dimension ( : ), allocatable :: &
        Solution_1D, &
        Y_1D
      type ( VariableGroupForm ), dimension ( :, : ), allocatable :: &
        K_1D
      type ( CurrentPointerForm ), dimension ( : ), allocatable :: &
        Current_1D
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_ASC_1D
    procedure, private, pass :: &
      Compute_C_ASC_1D
    generic, public :: &
      Compute => Compute_C_ASC_1D
    procedure, public, pass :: &
      FinalizeTemplate_C_ASC_1D
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


  subroutine InitializeTemplate_C_ASC_1D ( S, NameSuffix, A, B, C )

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

    call S % InitializeTemplate_C_ASC ( NameSuffix, A, B, C )

  end subroutine InitializeTemplate_C_ASC_1D


  subroutine Compute_C_ASC_1D &
               ( S, Current_ASC_1D, Time, TimeStep, UseLimiterParameterOption )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S
    class ( Current_ASC_ElementForm ), dimension ( : ), intent ( inout ), &
      target :: &
        Current_ASC_1D
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    logical ( KDL ), dimension ( : ), intent ( in ), optional :: &
      UseLimiterParameterOption

    integer ( KDI ) :: &
      iC  !-- iCurrent

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    S % nCurrents = size ( Current_ASC_1D )

    allocate ( S % UseLimiterParameter_C_1D ( S % nCurrents ) )
    S % UseLimiterParameter_C_1D = .true.
    if ( present ( UseLimiterParameterOption ) ) &
      S % UseLimiterParameter_C_1D = UseLimiterParameterOption

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
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_C_ASC_1D', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    call AllocateStorage ( S )
    call S % LoadSolution ( S % Solution_1D, S % Current_1D )

    ! call S % ComputeTemplate ( Time, TimeStep )

    call S % StoreSolution ( S % Current_1D, S % Solution_1D )
    call DeallocateStorage ( S )

    deallocate ( S % Current_1D )
    nullify ( S % Grid )
    deallocate ( S % UseLimiterParameter_C_1D )
    S % nCurrents = 0

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_C_ASC_1D


  impure elemental subroutine FinalizeTemplate_C_ASC_1D ( S )

    class ( Step_RK_C_ASC_1D_Template ), intent ( inout ) :: &
      S

    call S % FinalizeTemplate_C_ASC ( )

  end subroutine FinalizeTemplate_C_ASC_1D


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

      if ( allocated ( S % BoundaryFluence_CSL_C_1D ) ) &
        deallocate ( S % BoundaryFluence_CSL_C_1D )
      allocate ( S % BoundaryFluence_CSL_C_1D ( S % nCurrents ) )

      do iC = 1, S % nCurrents
        call S % AllocateBoundaryFluence &
               ( Grid, S % Current_1D ( iC ) % Pointer % N_CONSERVED, &
                 S % BoundaryFluence_CSL_C_1D ( iC ) % Array )
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
