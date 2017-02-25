!-- Step_RK_C_ASC is a template for a RungeKutta time step of a
!   conserved current.

module Step_RK_C_ASC__Template

  !-- Step_RungeKutta_Current_AtlasSingleChart_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Operations
  use Fields
  use Increments
  use Step_RK__Template

  implicit none
  private

!     type, private :: ApplyDivergencePointer
!       procedure ( ApplyDivergence ), pointer, nopass :: &
!         Pointer => ApplyDivergence
!     end type ApplyDivergencePointer

!     type, private :: ApplySourcesPointer
!       procedure ( AS ), pointer, nopass :: &
!         Pointer => null ( )
!     end type ApplySourcesPointer

!     type, private :: ApplyRelaxationPointer
!       procedure ( AR ), pointer, nopass :: &
!         Pointer => null ( )
!     end type ApplyRelaxationPointer

  type, public, extends ( Step_RK_Template ), abstract :: &
    Step_RK_C_ASC_Template
      integer ( KDI ) :: &
        iTimerGhost!, &
!         iGeometryValue
      type ( Real_1D_Form ), dimension ( : ), allocatable :: &
        dLogVolumeJacobian_dX
!       type ( Real_3D_2D_Form ), dimension ( : ), allocatable :: &
!         BoundaryFluence_CSL
      type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
        BoundaryFluence_CSL_C
      logical ( KDL ) :: &
        UseLimiterParameter_C
!       logical ( KDL ), dimension ( : ), allocatable :: &
!         UseLimiterParameter
      type ( VariableGroupForm ), allocatable :: &
        Solution_C, &
        Y_C
      type ( VariableGroupForm ), dimension ( : ), allocatable :: &
        K_C
!       class ( GeometryFlatForm ), pointer :: &
!         Geometry => null ( )
      class ( * ), pointer :: &
        Grid => null ( )
      class ( CurrentTemplate ), pointer :: &
        Current_C => null ( )
!       type ( CurrentPointerForm ), dimension ( : ), pointer :: &
!         Current_1D => null ( )
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
!       type ( ApplyDivergencePointer ), dimension ( : ), allocatable :: &
!         ApplyDivergence_1D
!       type ( ApplySourcesPointer ), dimension ( : ), allocatable :: &
!         ApplySources_1D
!       type ( ApplyRelaxationPointer ), dimension ( : ), allocatable :: &
!         ApplyRelaxation_1D
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_ASC
!     procedure, private, pass :: &
!       Compute_0D
!     procedure, private, pass :: &
!       Compute_1D
!     generic, public :: &
!       Compute => Compute_0D, Compute_1D
    procedure, private, pass :: &
      Compute_C_ASC
    generic, public :: &
      Compute => Compute_C_ASC
    procedure, public, pass :: &
      FinalizeTemplate_C_ASC
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
    procedure, private, nopass :: &
      LoadSolution_C
    generic, public :: &
      LoadSolution => LoadSolution_C
    procedure, private, pass ( S ) :: &
      StoreSolution_C
    generic, public :: &
      StoreSolution => StoreSolution_C
!     procedure, private, pass :: &
!       ComputeIncrement
!     procedure, public, pass :: &
!       SetDivergence
!     procedure, public, pass :: &
!       ClearDivergence
    procedure, public, pass :: &
      ComputeStage_C
    procedure, private, pass :: &
      AllocateBoundaryFluence_CSL
    generic, public :: &
      AllocateBoundaryFluence => AllocateBoundaryFluence_CSL
  end type Step_RK_C_ASC_Template

  abstract interface 

    subroutine AS ( S, Increment, Current, TimeStep )
      use Basics
      use Fields
      import Step_RK_C_ASC_Template
      class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
        S
      type ( VariableGroupForm ), intent ( inout ), target :: &
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
      import Step_RK_C_ASC_Template
      class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
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

    private :: &
      AllocateStorage, &
      DeallocateStorage

contains


  subroutine InitializeTemplate_C_ASC ( S, NameSuffix, A, B, C )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
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
      S % Type = 'a Step_RK_C_ASC' 

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

  end subroutine InitializeTemplate_C_ASC


!   subroutine Compute_0D &
!                ( S, Current, Grid, Time, TimeStep, GeometryOption, &
!                  UseLimiterParameterOption, iGeometryValueOption )

!     class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
!       S
!     class ( CurrentTemplate ), intent ( inout ), target :: &
!       Current
!     class ( * ), intent ( inout ), target :: &
!       Grid
!     real ( KDR ), intent ( in ) :: &
!       Time, &
!       TimeStep
!     class ( GeometryFlatForm ), intent ( in ), target, optional :: &
!       GeometryOption
!     logical ( KDL ), intent ( in ), optional :: &
!       UseLimiterParameterOption
!     integer ( KDI ), intent ( in ), optional :: &
!       iGeometryValueOption

!     logical ( KDL ), dimension ( 1 ) :: &
!       UseLimiterParameter
!     type ( CurrentPointerForm ), dimension ( 1 ) :: &
!       Current_1D

!     UseLimiterParameter = .true.
!     if ( present ( UseLimiterParameterOption ) ) &
!       UseLimiterParameter = UseLimiterParameterOption

!     Current_1D ( 1 ) % Pointer => Current

!     allocate ( S % ApplyDivergence_1D ( 1 ) )
!     allocate ( S % ApplySources_1D ( 1 ) )
!     allocate ( S % ApplyRelaxation_1D ( 1 ) )
!     S % ApplyDivergence_1D ( 1 ) % Pointer => S % ApplyDivergence
!     S % ApplySources_1D ( 1 )    % Pointer => S % ApplySources
!     S % ApplyRelaxation_1D ( 1 ) % Pointer => S % ApplyRelaxation

!     call S % Compute_1D &
!            ( Current_1D, Grid, Time, TimeStep, &
!              GeometryOption = GeometryOption, &
!              UseLimiterParameterOption = UseLimiterParameter, &
!              iGeometryValueOption = iGeometryValueOption )

!     deallocate ( S % ApplyRelaxation_1D )
!     deallocate ( S % ApplySources_1D )
!     deallocate ( S % ApplyDivergence_1D )

!   end subroutine Compute_0D


!   subroutine Compute_1D &
!                ( S, Current_1D, Grid, Time, TimeStep, GeometryOption, &
!                  UseLimiterParameterOption, iGeometryValueOption )

!     class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
!       S
!     type ( CurrentPointerForm ), dimension ( : ), intent ( inout ), &
!       target :: &
!         Current_1D
!     class ( * ), intent ( inout ), target :: &
!       Grid
!     real ( KDR ), intent ( in ) :: &
!       Time, &
!       TimeStep
!     class ( GeometryFlatForm ), intent ( in ), target, optional :: &
!       GeometryOption
!     logical ( KDL ), intent ( in ), dimension ( : ), optional :: &
!       UseLimiterParameterOption
!     integer ( KDI ), intent ( in ), optional :: &
!       iGeometryValueOption

!     integer ( KDI ) :: &
!       iF, &  !-- iField
!       iG     !-- iGroup
!     integer ( KDI ), dimension ( size ( Current_1D ) ) :: &
!       nValues, &
!       nEquations
!     type ( VariableGroupForm ), dimension ( size ( Current_1D ) ) :: &
!       Solution
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     associate &
!       ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
!     call Timer % Start ( )

!     S % Current_1D => Current_1D
!     S % Grid => Grid

!     associate ( nGroups => size ( Current_1D ) )

!     S % Geometry => null ( )
!     if ( present ( GeometryOption ) .and. present ( iGeometryValueOption ) ) &
!     then
!       S % Geometry => GeometryOption
!       S % iGeometryValue = iGeometryValueOption
!     end if

!     if ( .not. allocated ( S % ApplyDivergence_1D ) ) &
!       allocate ( S % ApplyDivergence_1D ( nGroups ) )
!     if ( .not. allocated ( S % ApplySources_1D ) ) &
!       allocate ( S % ApplySources_1D ( nGroups ) )
!     if ( .not. allocated ( S % ApplyRelaxation_1D ) ) &
!       allocate ( S % ApplyRelaxation_1D ( nGroups ) )

!     !-- Allocate Solution and initialize from Current_1D

!     do iG = 1, nGroups
!       associate &
!         ( C   => Current_1D ( iG ) % Pointer )
!       associate &
!         ( iaC => C % iaConserved, &
!           nE  => C % N_CONSERVED, &  !-- nEquations
!           nV  => C % nValues )

!       call Solution ( iG ) % Initialize ( [ nV, nE ] )
!       do iF = 1, C % N_CONSERVED
!         associate &
!           ( CV => C % Value ( :, iaC ( iF ) ), &
!             SV => Solution ( iG ) % Value ( :, iF ) )
!         call Copy ( CV, SV )
!         end associate !-- CV, etc.
!       end do !-- iF

!       end associate !-- iaC, etc.
!       end associate !-- C, etc.
!     end do !-- iG

!     !-- Compute Solution

!     call S % SetDivergence ( UseLimiterParameterOption )
!     call S % ComputeTemplate ( Solution, Time, TimeStep )
!     call S % ClearDivergence ( )

!     !-- Copy Solution to Current_1D

!     do iG = 1, nGroups
!       associate ( C => Current_1D ( iG ) % Pointer )
!       associate ( iaC => C % iaConserved )

!       do iF = 1, C % N_CONSERVED
!         associate &
!           ( SV => Solution ( iG ) % Value ( :, iF ), &
!             CV => C % Value ( :, iaC ( iF ) ) )
!         call Copy ( SV, CV )
!         end associate !-- YV, etc.
!       end do !-- iF

!       if ( associated ( S % Geometry ) ) then
!         call C % ComputeFromConserved ( S % iGeometryValue, S % Geometry )
!       else
!         select type ( Grid )
!         class is ( Chart_SL_Template )
!           G => Grid % Geometry ( )
!         class default
!           call Show ( 'Grid type not found', CONSOLE % ERROR )
!           call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
!           call Show ( 'Compute', 'subroutine', CONSOLE % ERROR ) 
!           call PROGRAM_HEADER % Abort ( )
!         end select !-- Grid
!         call C % ComputeFromConserved ( G )
!       end if

!       end associate !-- iaC
!       end associate !-- C
!     end do !-- iG

!     end associate !-- nGroups

!     nullify ( S % Grid )
!     nullify ( S % Current_1D )
!     nullify ( S % Geometry )

!     call Timer % Stop ( )
!     end associate !-- Timer

!   end subroutine Compute_1D


  subroutine Compute_C_ASC &
               ( S, Current_ASC, Time, TimeStep, UseLimiterParameterOption )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Current_ASC_Template ), intent ( inout ), target :: &
      Current_ASC
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameterOption

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    S % UseLimiterParameter_C = .true.
    if ( present ( UseLimiterParameterOption ) ) &
      S % UseLimiterParameter_C = UseLimiterParameterOption

    select type ( Chart => Current_ASC % Atlas_SC % Chart )
    class is ( Chart_SL_Template )
      S % Grid => Chart
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_C', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    S % Current_C => Current_ASC % Current ( )
    call AllocateStorage ( S, S % Current_C )
    call S % LoadSolution ( S % Solution_C, S % Current_C )

    call S % ComputeTemplate ( Time, TimeStep )

    call S % StoreSolution ( S % Current_C, S % Solution_C )
    call DeallocateStorage ( S )
    S % Current_C => null ( )
    S % Grid      => null ( )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_C_ASC


  impure elemental subroutine FinalizeTemplate_C_ASC ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % IncrementDamping ) ) &
      deallocate ( S % IncrementDamping )
    if ( allocated ( S % IncrementDivergence ) ) &
      deallocate ( S % IncrementDivergence )

    call S % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C_ASC


  subroutine InitializeIntermediate ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    associate &
      ( SV => S % Solution_C % Value, &
        YV => S % Y_C % Value )

    call Copy ( SV, YV )

    end associate !-- SV, etc.

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    associate &
      ( YV => S % Y_C % Value, &
        KV => S % K_C ( iK ) % Value )

    call MultiplyAdd ( YV, KV, A )

    end associate !-- YV, etc.

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, Time, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeIncrement ) )
    call Timer % Start ( )

    associate &
      ( C  => S % Current_C, &
        K  => S % K_C ( iStage ), &
        BF => S % BoundaryFluence_CSL_C, &
        Y  => S % Y_C )

    call S % ComputeStage_C &
           ( C, K, BF, Y, S % UseLimiterParameter_C, TimeStep, iStage )

    end associate !-- C, etc.

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, iS )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    associate &
      ( SV => S % Solution_C % Value, &
        KV => S % K_C ( iS ) % Value )

    call MultiplyAdd ( SV, KV, B )

    end associate !-- SV, etc.

  end subroutine IncrementSolution


  subroutine LoadSolution_C ( Solution_C, C )

    type ( VariableGroupForm ), intent ( inout ) :: &
      Solution_C
    class ( CurrentTemplate ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iF  !-- iField

    associate ( iaC => C % iaConserved )
    do iF = 1, C % N_CONSERVED
      associate &
        ( CV => C % Value ( :, iaC ( iF ) ), &
          SV => Solution_C % Value ( :, iF ) )
      call Copy ( CV, SV )
      end associate !-- CV, etc.
    end do !-- iF
    end associate !-- iaC

  end subroutine LoadSolution_C


  subroutine StoreSolution_C ( C, S, Solution_C )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( in ) :: &
      Solution_C

    integer ( KDI ) :: &
      iF  !-- iField
    class ( GeometryFlatForm ), pointer :: &
      G

    associate ( iaC => C % iaConserved )
    do iF = 1, C % N_CONSERVED
      associate &
        ( SV => Solution_C % Value ( :, iF ), &
          CV => C % Value ( :, iaC ( iF ) ) )
      call Copy ( SV, CV )
      end associate !-- YV, etc.
    end do !-- iF
    end associate !-- iaC

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )
      G => Grid % Geometry ( )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'StoreSolution_C', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid
    call C % ComputeFromConserved ( G )

    nullify ( G )
    
  end subroutine StoreSolution_C


!   subroutine ComputeIncrement ( S, K, Y, Time, TimeStep, iStage, iGroup )

!     class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
!       S
!     type ( VariableGroupForm ), dimension ( :, : ), intent ( inout ) :: &
!       K
!     type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
!       Y
!     real ( KDR ), intent ( in ) :: &
!       Time, &
!       TimeStep
!     integer ( KDI ), intent ( in ) :: &
!       iStage, &
!       iGroup

!     integer ( KDI ) :: &
!       iF  !-- iField
!     type ( VariableGroupForm ), allocatable :: &
!       DC  !-- DampingCoefficient
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     associate &
!       ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeIncrement ) )
!     call Timer % Start ( )

!     associate &
!       ( C   => S % Current_1D ( iGroup ) % Pointer, &
!         iaC => S % Current_1D ( iGroup ) % Pointer % iaConserved )

!     if ( iStage > 1 ) then

!       do iF = 1, C % N_CONSERVED
!         associate &
!           ( YV => Y ( iGroup ) % Value ( :, iF ), &
!             CV => C % Value ( :, iaC ( iF ) ) )
!         call Copy ( YV, CV )
!         end associate !-- YV, etc.
!       end do !-- iF

!       if ( associated ( S % Geometry ) ) then
!         call C % ComputeFromConserved ( S % iGeometryValue, S % Geometry )
!       else
!         select type ( Grid => S % Grid )
!         class is ( Chart_SL_Template )
!           G => Grid % Geometry ( )
!         class default
!           call Show ( 'Grid type not found', CONSOLE % ERROR )
!           call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
!           call Show ( 'ComputeIncrement', 'subroutine', CONSOLE % ERROR ) 
!           call PROGRAM_HEADER % Abort ( )
!         end select !-- Grid
!         call C % ComputeFromConserved ( G )
!       end if

!     end if !-- iStage > 1

!     associate ( KG => K ( iGroup, iStage ) )

!     call Clear ( KG % Value )

!     !-- Divergence
!     if ( associated ( S % ApplyDivergence_1D ( iGroup ) % Pointer ) ) &
!       call S % ApplyDivergence_1D ( iGroup ) % Pointer &
!              ( S, KG, C, TimeStep, iStage, iGroup )

!     !-- Other explicit sources
!     if ( associated ( S % ApplySources_1D ( iGroup ) % Pointer ) ) &
!       call S % ApplySources_1D ( iGroup ) % Pointer ( S, KG, C, TimeStep )

!     !-- Relaxation
!     if ( associated ( S % ApplyRelaxation_1D ( iGroup ) % Pointer ) ) then
!       associate ( ID => S % IncrementDamping )
!       allocate ( DC )
!       call DC % Initialize ( shape ( KG % Value ), ClearOption = .true. )
!       call S % ApplyRelaxation_1D ( iGroup ) % Pointer &
!              ( S, KG, DC, C, TimeStep )
!       call ID % Compute ( KG, C, KG, DC, TimeStep )
!       deallocate ( DC )
!       end associate !-- iD
!     end if

!     if ( associated ( S % ApplyDivergence ) ) then
!       select type ( Grid => S % Grid )
!       class is ( Chart_SLD_Form )
!         associate ( TimerGhost => PROGRAM_HEADER % Timer ( S % iTimerGhost ) )
!         call TimerGhost % Start ( )
!         call Grid % ExchangeGhostData ( KG )
!         call TimerGhost % Stop ( )
!         end associate !-- TimerGhost
!       end select !-- Grid
!     end if !-- ApplyDivergence

!     end associate !-- KG
!     end associate !-- C, etc.

!     call Timer % Stop ( )
!     end associate !-- Timer

!     nullify ( G )

!   end subroutine ComputeIncrement


!   subroutine SetDivergence ( S, UseLimiterParameterOption )

!     class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
!       S
!     logical ( KDL ), intent ( in ), dimension ( : ), optional :: &
!       UseLimiterParameterOption

!     integer ( KDI ) :: &
!       iD, jD, kD, &  !-- iDimension
!       iF, &  !-- iField
!       iG     !-- iGroup
!     integer ( KDI ), dimension ( 3 ) :: &
!       nSurface

!     associate ( nGroups => size ( S % Current_1D ) )

!     allocate ( S % UseLimiterParameter ( nGroups ) )
!     S % UseLimiterParameter = .true.
!     if ( present ( UseLimiterParameterOption ) ) &
!       S % UseLimiterParameter = UseLimiterParameterOption

!     select type ( Grid => S % Grid )
!     class is ( Chart_SLD_Form )

!       if ( allocated ( S % BoundaryFluence_CSL ) ) &
!         deallocate ( S % BoundaryFluence_CSL )
!       allocate ( S % BoundaryFluence_CSL ( nGroups ) )

!       do iG = 1, nGroups
!         associate &
!           ( C => Grid % Atlas % Connectivity, &
!             nDimensions => Grid % nDimensions, &
!             nConserved  => S % Current_1D ( iG ) % Pointer % N_CONSERVED, &
!             BF => S % BoundaryFluence_CSL ( iG ) )
!         call BF % Initialize ( [ nConserved, C % nFaces ] )
!         do iD = 1, nDimensions
!           jD = mod ( iD, 3 ) + 1
!           kD = mod ( jD, 3 ) + 1
!           nSurface ( iD ) = 1
!           nSurface ( jD ) = Grid % nCellsBrick ( jD ) 
!           nSurface ( kD ) = Grid % nCellsBrick ( kD )
!           do iF = 1, nConserved
!             call BF % Array ( iF, C % iaInner ( iD ) ) &
!                    % Initialize ( nSurface, ClearOption = .true. )
!             call BF % Array ( iF, C % iaOuter ( iD ) ) &
!                    % Initialize ( nSurface, ClearOption = .true. )
!           end do !-- iF
!         end do !-- iD
!         end associate !-- C, etc.

!       end do !-- iG

!       if ( ( trim ( Grid % CoordinateSystem ) == 'SPHERICAL' &
!              .or. trim ( Grid % CoordinateSystem ) == 'CYLINDRICAL' ) ) &
!       then

!         associate ( nValues => Grid % Geometry_CSL % nValues )
!         allocate ( S % dLogVolumeJacobian_dX ( 2 ) )
!         call S % dLogVolumeJacobian_dX ( 1 ) % Initialize ( nValues )
!         call S % dLogVolumeJacobian_dX ( 2 ) % Initialize ( nValues )
!         end associate !-- nValues

!         call S % IncrementDivergence % Set ( S % dLogVolumeJacobian_dX )

!       end if

!     class default
!       call Show ( 'Grid type not recognized', CONSOLE % ERROR )
!       call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
!       call Show ( 'SetDivergence', 'subroutine', CONSOLE % ERROR ) 
!       call PROGRAM_HEADER % Abort ( )
!     end select !-- Grid

!     end associate !-- nGroups

!   end subroutine SetDivergence


!   subroutine ClearDivergence ( S )

!     class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
!       S

!     if ( allocated ( S % dLogVolumeJacobian_dX ) ) &
!       deallocate ( S % dLogVolumeJacobian_dX )
        
!     deallocate ( S % UseLimiterParameter )

!   end subroutine ClearDivergence


  subroutine ComputeStage_C &
               ( S, C, K, BF, Y, UseLimiterParameter, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    type ( VariableGroupForm ), intent ( inout ) :: &
      K
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF
    type ( VariableGroupForm ), intent ( in ) :: &
      Y
    logical ( KDL ), intent ( in ) :: &
      UseLimiterParameter
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    type ( VariableGroupForm ), allocatable :: &
      DC  !-- DampingCoefficient

    if ( iStage > 1 ) &
      call S % StoreSolution ( C, Y )

    call Clear ( K % Value )

    !-- Divergence
    if ( associated ( S % ApplyDivergence ) ) &
      call S % ApplyDivergence &
             ( K, BF, C, UseLimiterParameter, TimeStep, iStage )

    !-- Other explicit sources
    if ( associated ( S % ApplySources ) ) &
      call S % ApplySources ( K, C, TimeStep )

    !-- Relaxation
    if ( associated ( S % ApplyRelaxation ) ) then
      associate ( ID => S % IncrementDamping )
      allocate ( DC )
      call DC % Initialize ( shape ( K % Value ), ClearOption = .true. )
      call S % ApplyRelaxation ( K, DC, C, TimeStep )
      call ID % Compute ( K, C, K, DC, TimeStep )
      deallocate ( DC )
      end associate !-- iD
    end if

    if ( associated ( S % ApplyDivergence ) ) then
      select type ( Grid => S % Grid )
      class is ( Chart_SLD_Form )
        associate ( TimerGhost => PROGRAM_HEADER % Timer ( S % iTimerGhost ) )
        call TimerGhost % Start ( )
        call Grid % ExchangeGhostData ( K )
        call TimerGhost % Stop ( )
        end associate !-- TimerGhost
      end select !-- Grid
    end if !-- ApplyDivergence

  end subroutine ComputeStage_C


  subroutine AllocateBoundaryFluence_CSL ( S, CSL, nEquations, BF )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      nEquations
    type ( Real_3D_Form ), dimension ( :, : ), intent ( out ), allocatable :: &
      BF

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension
      iE             !-- iEquation
    integer ( KDI ), dimension ( 3 ) :: &
      nSurface

      associate &
        ( C => CSL % Atlas % Connectivity, &
          nDimensions => CSL % nDimensions )

      if ( allocated ( BF ) ) &
        deallocate ( BF )
      allocate ( BF ( nEquations, C % nFaces ) )

      do iD = 1, nDimensions
        jD = mod ( iD, 3 ) + 1
        kD = mod ( jD, 3 ) + 1
        nSurface ( iD ) = 1
        select type ( CSL )
        class is ( Chart_SLL_Form )
          nSurface ( jD ) = CSL % nCells ( jD ) 
          nSurface ( kD ) = CSL % nCells ( kD )
        class is ( Chart_SLD_Form )
          nSurface ( jD ) = CSL % nCellsBrick ( jD ) 
          nSurface ( kD ) = CSL % nCellsBrick ( kD )
        end select !-- CSL
        do iE = 1, nEquations
          call BF ( iE, C % iaInner ( iD ) ) &
                 % Initialize ( nSurface, ClearOption = .true. )
          call BF ( iE, C % iaOuter ( iD ) ) &
                 % Initialize ( nSurface, ClearOption = .true. )
        end do !-- iE
      end do !-- iD
      end associate !-- C, etc.

  end subroutine AllocateBoundaryFluence_CSL


  subroutine ApplyDivergence &
               ( S, Increment, BoundaryFluence_CSL, Current, &
                 UseLimiterParameter, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BoundaryFluence_CSL
    class ( CurrentTemplate ), intent ( in ) :: &
      Current
    logical ( KDL ), intent ( in ) :: &
      UseLimiterParameter
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    associate ( ID => S % IncrementDivergence )
    call ID % Set ( UseLimiterParameter )
    call ID % Set ( BoundaryFluence_CSL )
    call ID % Set ( S % dLogVolumeJacobian_dX )
    call ID % Set ( Weight_RK = S % B ( iStage ) )
    call ID % Compute ( Increment, S % Grid, Current, TimeStep )
    call ID % Clear ( )
    end associate !-- ID

  end subroutine ApplyDivergence


  subroutine AllocateStorage ( S, C )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( CurrentTemplate ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iS  !-- iStage

    associate &
      ( nE  => C % N_CONSERVED, &  !-- nEquations
        nV  => C % nValues )

    allocate ( S % Solution_C )
    allocate ( S % Y_C )
    allocate ( S % K_C ( S % nStages ) )
    call S % Solution_C % Initialize ( [ nV, nE ] )
    call S % Y_C % Initialize ( [ nV, nE ] )
    do iS = 1, S % nStages
      call S % K_C ( iS ) % Initialize ( [ nV, nE ] )
    end do !-- iS

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )

      if ( ( trim ( Grid % CoordinateSystem ) == 'SPHERICAL' &
             .or. trim ( Grid % CoordinateSystem ) == 'CYLINDRICAL' ) ) &
      then
        allocate ( S % dLogVolumeJacobian_dX ( 2 ) )
        call S % dLogVolumeJacobian_dX ( 1 ) % Initialize ( nV )
        call S % dLogVolumeJacobian_dX ( 2 ) % Initialize ( nV )
      end if

      call S % AllocateBoundaryFluence ( Grid, nE, S % BoundaryFluence_CSL_C )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end associate !-- nE, etc.

  end subroutine AllocateStorage


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    !-- BoundaryFluence not deallocated here

    if ( allocated ( S % dLogVolumeJacobian_dX ) ) &
      deallocate ( S % dLogVolumeJacobian_dX )
    if ( allocated ( S % K_C ) ) &
      deallocate ( S % K_C )
    if ( allocated ( S % Y_C ) ) &
      deallocate ( S % Y_C )
    if ( allocated ( S % Solution_C ) ) &
      deallocate ( S % Solution_C )
      
  end subroutine DeallocateStorage


end module Step_RK_C_ASC__Template
