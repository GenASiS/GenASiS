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

  type, public :: ApplyDivergence_C_Pointer
    procedure ( ApplyDivergence_C ), pointer, nopass :: &
      Pointer => ApplyDivergence_C
  end type ApplyDivergence_C_Pointer

  type, public :: ApplySources_C_Pointer
    procedure ( AS ), pointer, nopass :: &
      Pointer => null ( )
  end type ApplySources_C_Pointer

  type, public :: ApplyRelaxation_C_Pointer
    procedure ( AR ), pointer, nopass :: &
      Pointer => null ( )
  end type ApplyRelaxation_C_Pointer

  type, public, extends ( Step_RK_Template ), abstract :: &
    Step_RK_C_ASC_Template
      integer ( KDI ) :: &
        iTimerGhost!, &
!         iGeometryValue
      type ( Real_1D_Form ), dimension ( : ), allocatable :: &
        dLogVolumeJacobian_dX
      type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
        BoundaryFluence_CSL
      type ( VariableGroupForm ), allocatable :: &
        Solution, &
        Y
      type ( VariableGroupForm ), dimension ( : ), allocatable :: &
        K
!       class ( GeometryFlatForm ), pointer :: &
!         Geometry => null ( )
      class ( ChartTemplate ), pointer :: &
        Chart => null ( )
      class ( CurrentTemplate ), pointer :: &
        Current => null ( )
      class ( Current_ASC_Template ), pointer :: &
        Current_ASC => null ( )
      type ( IncrementDivergence_FV_Form ), allocatable :: &
        IncrementDivergence_C
      type ( IncrementDampingForm ), allocatable :: &
        IncrementDamping
      procedure ( ApplyDivergence_C ), pointer, pass :: &
        ApplyDivergence_C => null ( )
      procedure ( AS ), pointer, pass :: &
        ApplySources_C => null ( ) 
      procedure ( AR ), pointer, pass :: &
        ApplyRelaxation_C => null ( )
      type ( ApplyDivergence_C_Pointer ) :: &
        ApplyDivergence
      type ( ApplySources_C_Pointer ) :: &
        ApplySources
      type ( ApplyRelaxation_C_Pointer ) :: &
        ApplyRelaxation
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_ASC
    procedure, public, pass :: &
      AllocateStorage
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
    procedure, public, pass :: &
      ComputeStage_C_ASC
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
      InitializeIntermediate_C
    procedure, public, pass :: &
      IncrementIntermediate_C
    procedure, public, pass :: &
      ComputeStage_C
    procedure, public, pass :: &
      IncrementSolution_C
    procedure, public, pass :: &
      Allocate_RK_C
    procedure, public, pass :: &
      Deallocate_RK_C
    procedure, private, nopass :: &
      AllocateBoundaryFluence_CSL
    generic, public :: &
      AllocateBoundaryFluence => AllocateBoundaryFluence_CSL
    procedure, private, nopass :: &
      ClearBoundaryFluence_CSL
    generic, public :: &
      ClearBoundaryFluence => ClearBoundaryFluence_CSL
    procedure, public, pass :: &
      DeallocateBoundaryFluence_CSL
    procedure, public, pass :: &
      AllocateMetricDerivatives
    procedure, public, pass :: &
      DeallocateMetricDerivatives
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
      DeallocateStorage, &
      ApplyDivergence_C

contains


  subroutine InitializeTemplate_C_ASC &
               ( S, Current_ASC, NameSuffix, A, B, C, UseLimiterOption )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Current_ASC_Template ), intent ( in ), target :: &
      Current_ASC
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_ASC' 

    call S % InitializeTemplate ( NameSuffix, A, B, C )

    select type ( Chart => Current_ASC % Atlas_SC % Chart )
    class is ( Chart_SL_Template )
      S % Chart => Chart
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C_ASC', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    S % Current_ASC => Current_ASC
    S % Current     => Current_ASC % Current ( )

    call S % AllocateStorage ( )

    !-- IncrementDivergence_C

    allocate ( S % IncrementDivergence_C )
    associate ( ID => S % IncrementDivergence_C )
    select type ( Chart => Current_ASC % Atlas_SC % Chart )
    class is ( Chart_SL_Template )
      call ID % Initialize &
             ( Current_ASC % Chart, UseLimiterOption, &
               BoundaryFluence_CSL_Option = S % BoundaryFluence_CSL, &
               dLogVolumeJacobian_dX_Option = S % dLogVolumeJacobian_dX )
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C_ASC', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart
    end associate !-- ID

    !-- IncrementDamping
    allocate ( S % IncrementDamping )
    associate ( ID => S % IncrementDamping )
    call ID % Initialize ( S % Name )
    end associate !-- ID

    call PROGRAM_HEADER % AddTimer &
           ( 'GhostIncrement', S % iTimerGhost )

  end subroutine InitializeTemplate_C_ASC


  subroutine AllocateStorage ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    character ( LDL ) :: &
      CoordinateSystem

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerAllocateStep ) )
    call Timer % Start ( )

    call DeallocateStorage ( S )

    call S % Allocate_RK_C ( )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

      call S % AllocateBoundaryFluence &
             ( Chart, S % Current % N_CONSERVED, S % BoundaryFluence_CSL )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    call S % AllocateMetricDerivatives &
           ( S % Chart % CoordinateSystem, S % Current % nValues )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine AllocateStorage


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


  subroutine Compute_C_ASC ( S, Time, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
    call Timer % Start ( )

    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      call S % ClearBoundaryFluence ( S % BoundaryFluence_CSL )

    call S % LoadSolution ( S % Solution, S % Current )
    call S % ComputeTemplate ( Time, TimeStep )
    call S % StoreSolution ( S % Current, S % Solution )

    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      call S % Current_ASC % AccumulateBoundaryTally ( S % BoundaryFluence_CSL )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine Compute_C_ASC


  impure elemental subroutine FinalizeTemplate_C_ASC ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    call DeallocateStorage ( S )

    if ( allocated ( S % IncrementDamping ) ) &
      deallocate ( S % IncrementDamping )
    if ( allocated ( S % IncrementDivergence_C ) ) &
      deallocate ( S % IncrementDivergence_C )

    nullify ( S % Current )
    nullify ( S % Current_ASC )

    call S % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C_ASC


  subroutine InitializeIntermediate ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    call S % InitializeIntermediate_C ( )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    call S % IncrementIntermediate_C  ( A, iK )

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
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStage ) )
    call Timer % Start ( )

    call S % ComputeStage_C_ASC ( TimeStep, iStage )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeStage


  subroutine ComputeStage_C_ASC ( S, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    associate &
      ( C     => S % Current, &
        Chart => S % Chart, &
        K     => S % K ( iStage ), &
        Y     => S % Y )

    if ( iStage > 1 ) &
      call S % StoreSolution ( C, Y )

    call Clear ( K % Value )

    S % ApplyDivergence_C => S % ApplyDivergence % Pointer
    S % ApplySources_C    => S % ApplySources    % Pointer
    S % ApplyRelaxation_C => S % ApplyRelaxation % Pointer

    call S % ComputeStage_C &
           ( S % IncrementDivergence_C, C, S % Chart, K, TimeStep, iStage )

    S % ApplyRelaxation_C => null ( )
    S % ApplySources_C    => null ( )
    S % ApplyDivergence_C => null ( )

    end associate !-- C, etc.

  end subroutine ComputeStage_C_ASC


  subroutine IncrementSolution ( S, B, iS )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    call S % IncrementSolution_C ( B, iS )

  end subroutine IncrementSolution


  subroutine LoadSolution_C ( Solution, Current )

    type ( VariableGroupForm ), intent ( inout ) :: &
      Solution
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    integer ( KDI ) :: &
      iF  !-- iField

    associate ( iaC => Current % iaConserved )
    do iF = 1, Current % N_CONSERVED
      associate &
        ( CV => Current % Value ( :, iaC ( iF ) ), &
          SV => Solution % Value ( :, iF ) )
      call Copy ( CV, SV )
      end associate !-- CV, etc.
    end do !-- iF
    end associate !-- iaC

  end subroutine LoadSolution_C


  subroutine StoreSolution_C ( Current, S, Solution )

    class ( CurrentTemplate ), intent ( inout ) :: &
      Current
    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( in ) :: &
      Solution

    integer ( KDI ) :: &
      iF  !-- iField
    class ( GeometryFlatForm ), pointer :: &
      G

    associate ( iaC => Current % iaConserved )
    do iF = 1, Current % N_CONSERVED
      associate &
        ( SV => Solution % Value ( :, iF ), &
          CV => Current % Value ( :, iaC ( iF ) ) )
      call Copy ( SV, CV )
      end associate !-- YV, etc.
    end do !-- iF
    end associate !-- iaC

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )
      G => Chart % Geometry ( )
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'StoreSolution_C', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid
    call Current % ComputeFromConserved ( G )

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


  subroutine InitializeIntermediate_C ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    associate &
      ( SV => S % Solution % Value, &
        YV => S % Y % Value )
    call Copy ( SV, YV )
    end associate !-- SV, etc.

  end subroutine InitializeIntermediate_C


  subroutine IncrementIntermediate_C ( S, A, iK )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    associate &
      ( YV => S % Y % Value, &
        KV => S % K ( iK ) % Value )
    call MultiplyAdd ( YV, KV, A )
    end associate !-- YV, etc.

  end subroutine IncrementIntermediate_C


  subroutine ComputeStage_C ( S, ID, C, Chart, K, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID    
    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( ChartTemplate ), intent ( inout ) :: &
      Chart
    type ( VariableGroupForm ), intent ( inout ) :: &
      K
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    type ( VariableGroupForm ), allocatable :: &
      DC  !-- DampingCoefficient

    !-- Divergence
    if ( associated ( S % ApplyDivergence_C ) ) &
      call S % ApplyDivergence_C ( ID, K, TimeStep, iStage )

    !-- Other explicit sources
    if ( associated ( S % ApplySources_C ) ) &
      call S % ApplySources_C ( K, C, TimeStep )

    !-- Relaxation
    if ( associated ( S % ApplyRelaxation_C ) ) then
      associate ( ID => S % IncrementDamping )
      allocate ( DC )
      call DC % Initialize ( shape ( K % Value ), ClearOption = .true. )
      call S % ApplyRelaxation_C ( K, DC, C, TimeStep )
      call ID % Compute ( K, C, K, DC, TimeStep )
      deallocate ( DC )
      end associate !-- iD
    end if

    if ( associated ( S % ApplyDivergence_C ) ) then
      select type ( Chart )
      class is ( Chart_SLD_Form )
        associate ( TimerGhost => PROGRAM_HEADER % Timer ( S % iTimerGhost ) )
        call TimerGhost % Start ( )
        call Chart % ExchangeGhostData ( K )
        call TimerGhost % Stop ( )
        end associate !-- TimerGhost
      end select !-- Grid
    end if !-- ApplyDivergence

  end subroutine ComputeStage_C


  subroutine IncrementSolution_C ( S, B, iS )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    associate &
      ( SV => S % Solution % Value, &
        KV => S % K ( iS ) % Value )
    call MultiplyAdd ( SV, KV, B )
    end associate !-- SV, etc.

  end subroutine IncrementSolution_C


  subroutine Allocate_RK_C ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iStage

    allocate ( S % Solution )
    allocate ( S % Y )
    allocate ( S % K ( S % nStages ) )

    associate &
      ( nEquations => S % Current % N_CONSERVED, &
        nValues    => S % Current % nValues )

    call S % Solution % Initialize ( [ nValues, nEquations ] )
    call S % Y % Initialize ( [ nValues, nEquations ] )
    do iS = 1, S % nStages
      call S % K ( iS ) % Initialize ( [ nValues, nEquations ] )
    end do !-- iS

    end associate !-- nEquations, etc.

  end subroutine Allocate_RK_C


  subroutine Deallocate_RK_C ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % K ) ) &
      deallocate ( S % K )
    if ( allocated ( S % Y ) ) &
      deallocate ( S % Y )
    if ( allocated ( S % Solution ) ) &
      deallocate ( S % Solution )
      
  end subroutine Deallocate_RK_C


  subroutine AllocateBoundaryFluence_CSL ( CSL, nEquations, BF )

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


  subroutine ClearBoundaryFluence_CSL ( BF )

    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF

    integer ( KDI ) :: &
      iE, &  !-- iEquation
      iC     !-- iConnectivity

    do iC = 1, size ( BF, dim = 2 )
      do iE = 1, size ( BF, dim = 1 )
        call Clear ( BF ( iE, iC ) % Value )
      end do !-- iE
    end do !-- iC

  end subroutine ClearBoundaryFluence_CSL


  subroutine DeallocateBoundaryFluence_CSL ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      deallocate ( S % BoundaryFluence_CSL )
    
  end subroutine DeallocateBoundaryFluence_CSL

  
  subroutine AllocateMetricDerivatives ( S, CoordinateSystem, nValues )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    integer ( KDI ), intent ( in ) :: &
      nValues

    if ( ( trim ( CoordinateSystem ) == 'SPHERICAL' &
           .or. trim ( CoordinateSystem ) == 'CYLINDRICAL' ) ) &
    then
      allocate ( S % dLogVolumeJacobian_dX ( 2 ) )
      call S % dLogVolumeJacobian_dX ( 1 ) % Initialize ( nValues )
      call S % dLogVolumeJacobian_dX ( 2 ) % Initialize ( nValues )
    end if

  end subroutine AllocateMetricDerivatives


  subroutine DeallocateMetricDerivatives ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    if ( allocated ( S % dLogVolumeJacobian_dX ) ) &
      deallocate ( S % dLogVolumeJacobian_dX )

  end subroutine DeallocateMetricDerivatives


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( S % iTimerAllocateStep ) )
    call Timer % Start ( )

    call S % DeallocateBoundaryFluence_CSL ( )
    call S % DeallocateMetricDerivatives ( )
    call S % Deallocate_RK_C ( )

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine DeallocateStorage


  subroutine ApplyDivergence_C ( S, ID, Increment, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call ID % Compute ( Increment, TimeStep, Weight_RK = S % B ( iStage ) )

  end subroutine ApplyDivergence_C


end module Step_RK_C_ASC__Template
