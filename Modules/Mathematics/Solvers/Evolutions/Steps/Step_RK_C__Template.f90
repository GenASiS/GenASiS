module Step_RK_C__Template

  !-- Step_RungeKutta_Current_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Fields
  use Increments
  use Step_RK__Template

  implicit none
  private

    type, private :: ApplyDivergencePointer
      procedure ( ApplyDivergence ), pointer, nopass :: &
        Pointer => ApplyDivergence
    end type ApplyDivergencePointer

    type, private :: ApplySourcesPointer
      procedure ( AS ), pointer, nopass :: &
        Pointer => null ( )
    end type ApplySourcesPointer

    type, private :: ApplyRelaxationPointer
      procedure ( AR ), pointer, nopass :: &
        Pointer => null ( )
    end type ApplyRelaxationPointer

  type, public, extends ( Step_RK_Template ), abstract :: Step_RK_C_Template
    integer ( KDI ) :: &
      iTimerGhost, &
      iGeometryValue
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      dLogVolumeJacobian_dX
    type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
      BoundaryFluence_CSL
    logical ( KDL ), dimension ( : ), allocatable :: &
      UseLimiterParameter
    class ( GeometryFlatForm ), pointer :: &
      Geometry
    class ( * ), pointer :: &
      Grid => null ( )
    type ( CurrentPointerForm ), dimension ( : ), pointer :: &
      Current_1D => null ( )
    type ( IncrementDivergence_FV_Form ), allocatable :: &
      IncrementDivergence
    type ( IncrementDampingForm ), allocatable :: &
      IncrementDamping
    procedure ( ApplyDivergence ), pointer, pass :: &
      ApplyDivergence => ApplyDivergence
    procedure ( AS ), pointer, pass :: &
      ApplySources => null ( ) 
    procedure ( AR ), pointer, pass :: &
      ApplyRelaxation => null ( )
    type ( ApplyDivergencePointer ), dimension ( : ), allocatable :: &
      ApplyDivergence_1D
    type ( ApplySourcesPointer ), dimension ( : ), allocatable :: &
      ApplySources_1D
    type ( ApplyRelaxationPointer ), dimension ( : ), allocatable :: &
      ApplyRelaxation_1D
  contains
    procedure, public, pass :: &
      InitializeTemplate_C
    procedure, private, pass :: &
      Compute_0D
    procedure, private, pass :: &
      Compute_1D
    generic, public :: &
      Compute => Compute_0D, Compute_1D
    procedure, public, pass :: &
      FinalizeTemplate_C
    procedure, private, pass :: &
      ComputeIncrement
    procedure, public, pass :: &
      SetDivergence
    procedure, public, pass :: &
      ClearDivergence
  end type Step_RK_C_Template

  abstract interface 

    subroutine AS ( S, Increment, Current, TimeStep )
      use Basics
      use Fields
      import Step_RK_C_Template
      class ( Step_RK_C_Template ), intent ( in ) :: &
        S
      type ( VariableGroupForm ), intent ( inout ) :: &
        Increment
      class ( CurrentTemplate ), intent ( in ) :: &
        Current
      real ( KDR ), intent ( in ) :: &
        TimeStep
    end subroutine AS

    subroutine AR ( S, IncrementExplicit, DampingCoefficient, Current, &
                    TimeStep )
      use Basics
      use Fields
      import Step_RK_C_Template
      class ( Step_RK_C_Template ), intent ( in ) :: &
        S
      type ( VariableGroupForm ), intent ( inout ) :: &
        IncrementExplicit, &
        DampingCoefficient
      class ( CurrentTemplate ), intent ( in ) :: &
        Current
      real ( KDR ), intent ( in ) :: &
        TimeStep
    end subroutine AR
    
  end interface

contains


  subroutine InitializeTemplate_C ( S, NameSuffix, A, B, C )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
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
      S % Type = 'a Step_RK_C' 

    call S % InitializeTemplate ( NameSuffix, A, B, C )

    allocate ( S % IncrementDivergence )
    associate ( ID => S % IncrementDivergence )
    call ID % Initialize ( S % Name )
    end associate !-- ID

    allocate ( S % IncrementDamping )
    associate ( ID => S % IncrementDamping )
    call ID % Initialize ( S % Name )
    end associate !-- ID

    call PROGRAM_HEADER % AddTimer &
           ( 'GhostIncrement', S % iTimerGhost )

  end subroutine InitializeTemplate_C


  subroutine Compute_0D &
               ( S, Current, Grid, Time, TimeStep, GeometryOption, &
                 UseLimiterParameterOption, iGeometryValueOption )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    class ( CurrentTemplate ), intent ( inout ), target :: &
      Current
    class ( * ), intent ( inout ), target :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    class ( GeometryFlatForm ), intent ( in ), target, optional :: &
      GeometryOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameterOption
    integer ( KDI ), intent ( in ), optional :: &
      iGeometryValueOption

    logical ( KDL ), dimension ( 1 ) :: &
      UseLimiterParameter
    type ( CurrentPointerForm ), dimension ( 1 ) :: &
      Current_1D

    UseLimiterParameter = .true.
    if ( present ( UseLimiterParameterOption ) ) &
      UseLimiterParameter = UseLimiterParameterOption

    Current_1D ( 1 ) % Pointer => Current

    allocate ( S % ApplyDivergence_1D ( 1 ) )
    allocate ( S % ApplySources_1D ( 1 ) )
    allocate ( S % ApplyRelaxation_1D ( 1 ) )
    S % ApplyDivergence_1D ( 1 ) % Pointer => S % ApplyDivergence
    S % ApplySources_1D ( 1 )    % Pointer => S % ApplySources
    S % ApplyRelaxation_1D ( 1 ) % Pointer => S % ApplyRelaxation

    call S % Compute_1D &
           ( Current_1D, Grid, Time, TimeStep, &
             GeometryOption = GeometryOption, &
             UseLimiterParameterOption = UseLimiterParameter, &
             iGeometryValueOption = iGeometryValueOption )

    deallocate ( S % ApplyRelaxation_1D )
    deallocate ( S % ApplySources_1D )
    deallocate ( S % ApplyDivergence_1D )

  end subroutine Compute_0D


  subroutine Compute_1D &
               ( S, Current_1D, Grid, Time, TimeStep, GeometryOption, &
                 UseLimiterParameterOption, iGeometryValueOption )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    type ( CurrentPointerForm ), dimension ( : ), intent ( inout ), &
      target :: &
        Current_1D
    class ( * ), intent ( inout ), target :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    class ( GeometryFlatForm ), intent ( in ), target, optional :: &
      GeometryOption
    logical ( KDL ), intent ( in ), dimension ( : ), optional :: &
      UseLimiterParameterOption
    integer ( KDI ), intent ( in ), optional :: &
      iGeometryValueOption

    integer ( KDI ) :: &
      iF, &  !-- iField
      iG     !-- iGroup
    integer ( KDI ), dimension ( size ( Current_1D ) ) :: &
      nValues, &
      nEquations
    type ( VariableGroupForm ), dimension ( size ( Current_1D ) ) :: &
      Solution
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    S % Current_1D => Current_1D
    S % Grid => Grid

    associate ( nGroups => size ( Current_1D ) )

    if ( present ( GeometryOption ) .and. present ( iGeometryValueOption ) ) &
    then
      S % Geometry => GeometryOption
      S % iGeometryValue = iGeometryValueOption
    end if

    if ( .not. allocated ( S % ApplyDivergence_1D ) ) &
      allocate ( S % ApplyDivergence_1D ( nGroups ) )
    if ( .not. allocated ( S % ApplySources_1D ) ) &
      allocate ( S % ApplySources_1D ( nGroups ) )
    if ( .not. allocated ( S % ApplyRelaxation_1D ) ) &
      allocate ( S % ApplyRelaxation_1D ( nGroups ) )

    !-- Allocate Solution and initialize from Current_1D

    do iG = 1, nGroups
      associate &
        ( C   => Current_1D ( iG ) % Pointer )
      associate &
        ( iaC => C % iaConserved, &
          nE  => C % N_CONSERVED, &  !-- nEquations
          nV  => C % nValues )

      call Solution ( iG ) % Initialize ( [ nV, nE ] )
      do iF = 1, C % N_CONSERVED
        associate &
          ( CV => C % Value ( :, iaC ( iF ) ), &
            SV => Solution ( iG ) % Value ( :, iF ) )
        call Copy ( CV, SV )
        end associate !-- CV, etc.
      end do !-- iF

      end associate !-- iaC, etc.
      end associate !-- C, etc.
    end do !-- iG

    !-- Compute Solution

    call S % SetDivergence ( UseLimiterParameterOption )
    call S % ComputeTemplate ( Solution, Time, TimeStep )
    call S % ClearDivergence ( )

    !-- Copy Solution to Current_1D

    do iG = 1, nGroups
      associate ( C => Current_1D ( iG ) % Pointer )
      associate ( iaC => C % iaConserved )

      do iF = 1, C % N_CONSERVED
        associate &
          ( SV => Solution ( iG ) % Value ( :, iF ), &
            CV => C % Value ( :, iaC ( iF ) ) )
        call Copy ( SV, CV )
        end associate !-- YV, etc.
      end do !-- iF

      if ( associated ( S % Geometry ) ) then
        call C % ComputeFromConserved ( S % iGeometryValue, S % Geometry )
      else
        select type ( Grid )
        class is ( Chart_SL_Template )
          G => Grid % Geometry ( )
        class default
          call Show ( 'Grid type not found', CONSOLE % ERROR )
          call Show ( 'Step_RK_C__Form', 'module', CONSOLE % ERROR )
          call Show ( 'Compute', 'subroutine', CONSOLE % ERROR ) 
          call PROGRAM_HEADER % Abort ( )
        end select !-- Grid
        call C % ComputeFromConserved ( G )
      end if

      end associate !-- iaC
      end associate !-- C
    end do !-- iG

    end associate !-- nGroups

    nullify ( S % Grid )
    nullify ( S % Current_1D )
    nullify ( S % Geometry )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_1D


  impure elemental subroutine FinalizeTemplate_C ( S )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % IncrementDivergence ) ) &
      deallocate ( S % IncrementDivergence )
    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      deallocate ( S % BoundaryFluence_CSL )

    nullify ( S % Current_1D )
    nullify ( S % Grid )

    call S % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C


  subroutine ComputeIncrement ( S, K, Y, Time, TimeStep, iStage, iGroup )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    type ( VariableGroupForm ), dimension ( :, : ), intent ( inout ) :: &
      K
    type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      Y
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage, &
      iGroup

    integer ( KDI ) :: &
      iF  !-- iField
    type ( VariableGroupForm ), allocatable :: &
      DC  !-- DampingCoefficient
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeIncrement ) )
    call Timer % Start ( )

    associate &
      ( C   => S % Current_1D ( iGroup ) % Pointer, &
        iaC => S % Current_1D ( iGroup ) % Pointer % iaConserved )

    if ( iStage > 1 ) then

      do iF = 1, C % N_CONSERVED
        associate &
          ( YV => Y ( iGroup ) % Value ( :, iF ), &
            CV => C % Value ( :, iaC ( iF ) ) )
        call Copy ( YV, CV )
        end associate !-- YV, etc.
      end do !-- iF

      if ( associated ( S % Geometry ) ) then
        call C % ComputeFromConserved ( S % iGeometryValue, S % Geometry )
      else
        select type ( Grid => S % Grid )
        class is ( Chart_SL_Template )
          G => Grid % Geometry ( )
        class default
          call Show ( 'Grid type not found', CONSOLE % ERROR )
          call Show ( 'Step_RK_C__Form', 'module', CONSOLE % ERROR )
          call Show ( 'ComputeIncrement', 'subroutine', CONSOLE % ERROR ) 
          call PROGRAM_HEADER % Abort ( )
        end select !-- Grid
        call C % ComputeFromConserved ( G )
      end if

    end if !-- iStage > 1

    associate ( KG => K ( iGroup, iStage ) )

    call Clear ( KG % Value )

    !-- Divergence
    if ( associated ( S % ApplyDivergence_1D ( iGroup ) % Pointer ) ) &
      call S % ApplyDivergence_1D ( iGroup ) % Pointer &
             ( S, KG, C, TimeStep, iStage, iGroup )

    !-- Other explicit sources
    if ( associated ( S % ApplySources_1D ( iGroup ) % Pointer ) ) &
      call S % ApplySources_1D ( iGroup ) % Pointer ( S, KG, C, TimeStep )

    !-- Relaxation
    if ( associated ( S % ApplyRelaxation_1D ( iGroup ) % Pointer ) ) then
      associate ( ID => S % IncrementDamping )
      allocate ( DC )
      call DC % Initialize ( shape ( KG % Value ), ClearOption = .true. )
      call S % ApplyRelaxation_1D ( iGroup ) % Pointer &
             ( S, KG, DC, C, TimeStep )
      call ID % Compute ( KG, C, KG, DC, TimeStep )
      deallocate ( DC )
      end associate !-- iD
    end if

    if ( associated ( S % ApplyDivergence ) ) then
      select type ( Grid => S % Grid )
      class is ( Chart_SLD_Form )
        associate ( TimerGhost => PROGRAM_HEADER % Timer ( S % iTimerGhost ) )
        call TimerGhost % Start ( )
        call Grid % ExchangeGhostData ( KG )
        call TimerGhost % Stop ( )
        end associate !-- TimerGhost
      end select !-- Grid
    end if !-- ApplyDivergence

    end associate !-- KG
    end associate !-- C, etc.

    call Timer % Stop ( )
    end associate !-- Timer

    nullify ( G )

  end subroutine ComputeIncrement


  subroutine SetDivergence ( S, UseLimiterParameterOption )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    logical ( KDL ), intent ( in ), dimension ( : ), optional :: &
      UseLimiterParameterOption

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension
      iF, &  !-- iField
      iG     !-- iGroup
    integer ( KDI ), dimension ( 3 ) :: &
      nSurface

    associate ( nGroups => size ( S % Current_1D ) )

    allocate ( S % UseLimiterParameter ( nGroups ) )
    S % UseLimiterParameter = .true.
    if ( present ( UseLimiterParameterOption ) ) &
      S % UseLimiterParameter = UseLimiterParameterOption

    select type ( Grid => S % Grid )
    class is ( Chart_SLD_Form )

      if ( allocated ( S % BoundaryFluence_CSL ) ) &
        deallocate ( S % BoundaryFluence_CSL )
      allocate ( S % BoundaryFluence_CSL ( nGroups ) )

      do iG = 1, nGroups
        associate &
          ( C => Grid % Atlas % Connectivity, &
            nDimensions => Grid % nDimensions, &
            nConserved  => S % Current_1D ( iG ) % Pointer % N_CONSERVED, &
            BF => S % BoundaryFluence_CSL ( iG ) )
        call BF % Initialize ( [ nConserved, C % nFaces ] )
        do iD = 1, nDimensions
          jD = mod ( iD, 3 ) + 1
          kD = mod ( jD, 3 ) + 1
          nSurface ( iD ) = 1
          nSurface ( jD ) = Grid % nCellsBrick ( jD ) 
          nSurface ( kD ) = Grid % nCellsBrick ( kD )
          do iF = 1, nConserved
            call BF % Array ( iF, C % iaInner ( iD ) ) &
                   % Initialize ( nSurface, ClearOption = .true. )
            call BF % Array ( iF, C % iaOuter ( iD ) ) &
                   % Initialize ( nSurface, ClearOption = .true. )
          end do !-- iF
        end do !-- iD
        end associate !-- C, etc.

      end do !-- iG

      if ( ( trim ( Grid % CoordinateSystem ) == 'SPHERICAL' &
             .or. trim ( Grid % CoordinateSystem ) == 'CYLINDRICAL' ) ) &
      then

        associate ( nValues => Grid % Geometry_CSL % nValues )
        allocate ( S % dLogVolumeJacobian_dX ( 2 ) )
        call S % dLogVolumeJacobian_dX ( 1 ) % Initialize ( nValues )
        call S % dLogVolumeJacobian_dX ( 2 ) % Initialize ( nValues )
        end associate !-- nValues

        call S % IncrementDivergence % Set ( S % dLogVolumeJacobian_dX )

      end if

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C__Template', 'module', CONSOLE % ERROR )
      call Show ( 'SetDivergence', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end associate !-- nGroups

  end subroutine SetDivergence


  subroutine ClearDivergence ( S )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % dLogVolumeJacobian_dX ) ) &
      deallocate ( S % dLogVolumeJacobian_dX )
        
    deallocate ( S % UseLimiterParameter )

  end subroutine ClearDivergence


  subroutine ApplyDivergence &
               ( S, Increment, Current, TimeStep, iStage, iGroup )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Current
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage, &
      iGroup

    associate ( ID => S % IncrementDivergence )
    call ID % Set ( S % UseLimiterParameter ( iGroup ) )
    call ID % Set ( S % BoundaryFluence_CSL ( iGroup ) % Array )
    call ID % Set ( Weight_RK = S % B ( iStage ) )
    call ID % Compute ( Increment, S % Grid, Current, TimeStep )
    call ID % Clear ( )
    end associate !-- ID

  end subroutine ApplyDivergence


end module Step_RK_C__Template
