!-- Step_RK_C is a template for a RungeKutta time step of a conserved current.

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

  type, public, extends ( Step_RK_Template ), abstract :: Step_RK_C_Template
    integer ( KDI ) :: &
      iTimerGhost
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      dLogVolumeJacobian_dX
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      BoundaryFluence_CSL
    class ( * ), pointer :: &
      Grid => null ( )
    class ( CurrentTemplate ), pointer :: &
      Current => null ( )
    type ( IncrementDivergence_FV_Form ), allocatable :: &
      IncrementDivergence
    procedure ( AS ), pointer, pass :: &
      ApplySources => null ( ) 
    procedure ( AR ), pointer, pass :: &
      ApplyRelaxation => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate_C
    procedure, public, pass :: &
      Compute
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

    call PROGRAM_HEADER % AddTimer &
           ( 'GhostIncrement', S % iTimerGhost )

  end subroutine InitializeTemplate_C


  subroutine Compute ( S, Current, Grid, Time, TimeStep )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    class ( CurrentTemplate ), intent ( inout ), target :: &
      Current
    class ( * ), intent ( inout ), target :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    integer ( KDI ) :: &
      iF  !-- iField
    type ( VariableGroupForm ) :: &
      Solution
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    S % Current => Current
    S % Grid => Grid

    associate &
      ( C   => Current, &
        iaC => Current % iaConserved, &
        nE  => Current % N_CONSERVED, &  !-- nEquations
        nV  => Current % nValues )

    select type ( Grid )
    class is ( Chart_SL_Template )
      G => Grid % Geometry ( )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    call Solution % Initialize ( [ nV, nE ] )
    do iF = 1, C % N_CONSERVED
      associate &
        ( CV => C % Value ( :, iaC ( iF ) ), &
          SV => Solution % Value ( :, iF ) )
      call Copy ( CV, SV )
      end associate !-- CV, etc.
    end do !-- iF

    call S % SetDivergence ( )
    call S % ComputeTemplate ( Solution, Time, TimeStep )

    do iF = 1, C % N_CONSERVED
      associate &
        ( SV => Solution % Value ( :, iF ), &
          CV => C % Value ( :, iaC ( iF ) ) )
      call Copy ( SV, CV )
      end associate !-- YV, etc.
    end do !-- iF
    call C % ComputeFromConserved ( G )

    call S % ClearDivergence ( )

    end associate !-- C, etc.

    nullify ( S % Grid )
    nullify ( S % Current )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute


  impure elemental subroutine FinalizeTemplate_C ( S )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % IncrementDivergence ) ) &
      deallocate ( S % IncrementDivergence )
    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      deallocate ( S % BoundaryFluence_CSL )

    nullify ( S % Current )
    nullify ( S % Grid )

    call S % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C


  subroutine ComputeIncrement ( S, K, Y, Time, TimeStep, iStage )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ) :: &
      K
    type ( VariableGroupForm ), intent ( in ) :: &
      Y
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iF  !-- iField
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeIncrement ) )
    call Timer % Start ( )

    associate &
      ( C   => S % Current, &
        iaC => S % Current % iaConserved )

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )
      G => Grid % Geometry ( )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeIncrement', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    if ( iStage > 1 ) then
      do iF = 1, C % N_CONSERVED
        associate &
          ( YV => Y % Value ( :, iF ), &
            CV => C % Value ( :, iaC ( iF ) ) )
        call Copy ( YV, CV )
        end associate !-- YV, etc.
      end do !-- iF
      call C % ComputeFromConserved ( G )
    end if

    call Clear ( K % Value )

    !-- Divergence
    associate ( ID => S % IncrementDivergence )
    call ID % Set ( Weight_RK = S % B ( iStage ) )
    call ID % Compute ( K, S % Grid, C, TimeStep )
    end associate !-- ID

    !-- Other explicit sources
    if ( associated ( S % ApplySources ) ) &
      call S % ApplySources ( K, C, TimeStep )

    select type ( Grid => S % Grid )
    class is ( Chart_SLD_Form )
      associate ( TimerGhost => PROGRAM_HEADER % Timer ( S % iTimerGhost ) )
      call TimerGhost % Start ( )
      call Grid % ExchangeGhostData ( K )
      call TimerGhost % Stop ( )
      end associate !-- TimerGhost
    end select !-- Grid

    end associate !-- C, etc.

    call Timer % Stop ( )
    end associate !-- Timer

    nullify ( G )

  end subroutine ComputeIncrement


  subroutine SetDivergence ( S )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension
      iF      !-- iField
    integer ( KDI ), dimension ( 3 ) :: &
      nSurface

    select type ( Grid => S % Grid)
    class is ( Chart_SLD_Form )

      if ( allocated ( S % BoundaryFluence_CSL ) ) &
        deallocate ( S % BoundaryFluence_CSL )

      associate &
        ( C => Grid % Atlas % Connectivity, &
          nDimensions => Grid % nDimensions, &
          nConserved => S % Current % N_CONSERVED )
      allocate &
        ( S % BoundaryFluence_CSL ( nConserved, C % nFaces ) )
      do iD = 1, nDimensions
        jD = mod ( iD, 3 ) + 1
        kD = mod ( jD, 3 ) + 1
        nSurface ( iD ) = 1
        nSurface ( jD ) = Grid % nCellsBrick ( jD ) 
        nSurface ( kD ) = Grid % nCellsBrick ( kD )
        do iF = 1, S % Current % N_CONSERVED
          call S % BoundaryFluence_CSL ( iF, C % iaInner ( iD ) ) &
                 % Initialize ( nSurface, ClearOption = .true. )
          call S % BoundaryFluence_CSL ( iF, C % iaOuter ( iD ) ) &
                 % Initialize ( nSurface, ClearOption = .true. )
        end do !-- iF
      end do !-- iD
      end associate !-- C, etc.

      call S % IncrementDivergence % Set ( S % BoundaryFluence_CSL )

      if ( trim ( Grid % CoordinateSystem ) == 'SPHERICAL' &
           .or. trim ( Grid % CoordinateSystem ) == 'CYLINDRICAL' ) then

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

  end subroutine SetDivergence


  subroutine ClearDivergence ( S )

    class ( Step_RK_C_Template ), intent ( inout ) :: &
      S

    call S % IncrementDivergence % Clear ( )

    if ( allocated ( S % dLogVolumeJacobian_dX ) ) &
      deallocate ( S % dLogVolumeJacobian_dX )
        
  end subroutine ClearDivergence


end module Step_RK_C__Template
