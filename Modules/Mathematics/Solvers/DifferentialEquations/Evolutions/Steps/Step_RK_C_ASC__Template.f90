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
  use EvolutionBasics
  use Increments
  use Step_RK__Template

  implicit none
  private

  type, public :: ComputeConstraints_C_Pointer
    procedure ( CC ), pointer, nopass :: &
      Pointer => null ( )
  end type ComputeConstraints_C_Pointer

  type, public :: ApplyDivergence_C_Pointer
    procedure ( ApplyDivergence_C ), pointer, nopass :: &
      Pointer => null ( )
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
        iTimerStoreIntermediate = 0, &
        iTimerClearIncrement    = 0, &
        iTimerConstraints       = 0, &
        iTimerDivergence        = 0, &
        iTimerSources           = 0, &
        iTimerRelaxation        = 0, &
        iTimerGhost             = 0, &
        iTimerBoundaryFluence   = 0
!         iStrgeometryValue
      type ( Real_1D_Form ), dimension ( : ), allocatable :: &
        dLogVolumeJacobian_dX
      type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
        BoundaryFluence_CSL
      logical ( KDL ) :: &
        Allocated
      type ( StorageForm ), allocatable :: &
        Solution, &
        Y
      type ( StorageForm ), dimension ( : ), allocatable :: &
        K
!       class ( GeometryFlatForm ), pointer :: &
!         Geometry => null ( )
      class ( ChartTemplate ), pointer :: &
        Chart => null ( )
      class ( CurrentTemplate ), pointer :: &
        Current => null ( )
      class ( Current_ASC_Template ), pointer :: &
        Current_ASC => null ( )
      type ( StorageDivergenceForm ), allocatable :: &
        StorageDivergence_C
      type ( IncrementDivergence_FV_Form ), allocatable :: &
        IncrementDivergence_C
!      type ( IncrementDampingForm ), allocatable :: &
!        IncrementDamping
      procedure ( ApplyDivergence_C ), pointer, pass :: &
        ApplyDivergence_C => null ( )
      procedure ( AS ), pointer, pass :: &
        ApplySources_C => null ( ) 
      procedure ( AR ), pointer, pass :: &
        ApplyRelaxation_C => null ( )
      type ( ComputeConstraints_C_Pointer ) :: &
        ComputeConstraints
      type ( ApplyDivergence_C_Pointer ) :: &
        ApplyDivergence
      type ( ApplySources_C_Pointer ) :: &
        ApplySources
      type ( ApplyRelaxation_C_Pointer ) :: &
        ApplyRelaxation
      procedure ( CS ), pointer, pass :: &
        CoarsenSingularities => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_ASC
    procedure, public, pass :: &
      InitializeTimers
    procedure, public, pass :: &
      DeallocateStorage
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      FinalizeTemplate_C_ASC
    procedure, public, pass :: &
      SetUseLimiter
    procedure, public, pass :: &
      InitializeTimersStage
    procedure, private, pass :: &
      InitializeTimersDivergence
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
    procedure, public, pass :: &
      ComputeStage_C_ASC
    procedure, private, pass :: &
      IncrementSolution
    procedure, public, pass ( S ) :: &
      LoadSolution_C
    procedure, public, pass ( S ) :: &
      StoreSolution_C
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
    procedure, public, nopass :: &
      AllocateStorageDivergence
    procedure, public, nopass :: &
      DeallocateStorageDivergence
    procedure, private, nopass :: &
      AllocateBoundaryFluence_CSL
    generic, public :: &
      AllocateBoundaryFluence => AllocateBoundaryFluence_CSL
    procedure, private, nopass :: &
      ClearBoundaryFluence_CSL
    generic, public :: &
      ClearBoundaryFluence => ClearBoundaryFluence_CSL
    procedure, public, nopass :: &
      DeallocateBoundaryFluence_CSL
    procedure, public, pass :: &
      AllocateMetricDerivatives
    procedure, public, pass :: &
      DeallocateMetricDerivatives
  end type Step_RK_C_ASC_Template

  abstract interface 

    subroutine CC ( S )
      import Step_RK_C_ASC_Template
      class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
        S
    end subroutine CC

    subroutine AS ( S, Sources, Increment, Current, TimeStep, iStage )
      use Basics
      use Fields
      import Step_RK_C_ASC_Template
      class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
        S
      class ( Sources_C_Form ), intent ( inout ) :: &
        Sources
      type ( StorageForm ), intent ( inout ), target :: &
        Increment
      class ( CurrentTemplate ), intent ( in ) :: &
        Current
      real ( KDR ), intent ( in ) :: &
        TimeStep
      integer ( KDI ), intent ( in ) :: &
        iStage
    end subroutine AS

    subroutine AR ( S, Current, Sources, Increment, Chart, TimeStep, iStage, &
                    GeometryOption, iStrgeometryValueOption )
      use Basics
      use Manifolds
      use Fields
      import Step_RK_C_ASC_Template
      class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
        S
      class ( CurrentTemplate ), intent ( inout ) :: &
        Current
      class ( Sources_C_Form ), intent ( inout ) :: &
        Sources
      type ( StorageForm ), intent ( inout ) :: &
        Increment
      class ( ChartTemplate ), intent ( in ) :: &
        Chart
      real ( KDR ), intent ( in ) :: &
        TimeStep
      integer ( KDI ), intent ( in ) :: &
        iStage
      class ( GeometryFlatForm ), intent ( in ), optional :: &
        GeometryOption
      integer ( KDI ), intent ( in ), optional :: &
        iStrgeometryValueOption
    end subroutine AR
    
    subroutine CS ( S, Increment )
      use Basics
      import Step_RK_C_ASC_Template
      class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
        S
      class ( StorageForm ), intent ( inout ) :: &
        Increment
    end subroutine CS

  end interface

  public :: &
    ApplyDivergence_C

    private :: &
      AllocateStorage, &
      RecordDivergence

contains


  subroutine InitializeTemplate_C_ASC &
               ( S, I, Current_ASC, NameSuffix, A, B, C )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( IntegratorHeaderForm ), intent ( in ) :: &
      I
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

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_ASC' 

    call S % InitializeTemplate ( I, NameSuffix, A, B, C )

    select type ( Chart => Current_ASC % Atlas_SC % Chart )
    class is ( Chart_SL_Template )
      S % Chart => Chart
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C_ASC', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    S % Current_ASC => Current_ASC
    S % Current     => Current_ASC % Current ( )

    S % Allocated = .false.

    allocate ( S % StorageDivergence_C )
    allocate ( S % IncrementDivergence_C )
    associate ( ID => S % IncrementDivergence_C )
    call ID % Initialize ( S % Current_ASC % Chart )
    call ID % SetStorage ( S % StorageDivergence_C )
    end associate !-- ID

!    allocate ( S % IncrementDamping )
!    associate ( ID => S % IncrementDamping )
!    call ID % Initialize ( S % Name )
!    end associate !-- ID

    S % ApplyDivergence % Pointer => ApplyDivergence_C

  end subroutine InitializeTemplate_C_ASC


  subroutine InitializeTimers ( S, BaseLevel )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( S % iTimerStep > 0  &
         .or.  BaseLevel > PROGRAM_HEADER % TimerLevel ) &
      return

    call PROGRAM_HEADER % AddTimer &
           ( S % Name, S % iTimerStep, &
             Level = BaseLevel )
      call S % InitializeTimersTemplate ( BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'BoundaryFluence', S % iTimerBoundaryFluence, &
               Level = BaseLevel + 1 )

  end subroutine InitializeTimers


  subroutine DeallocateStorage ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    S % Allocated = .false.

    if ( .not. allocated ( S % IncrementDivergence_C ) ) &
      return

    associate &
      ( ID => S % IncrementDivergence_C, &
        SD => S % StorageDivergence_C )
    call S % DeallocateMetricDerivatives ( ID )
    call S % DeallocateBoundaryFluence_CSL ( ID, S % BoundaryFluence_CSL )
    call S % DeallocateStorageDivergence ( SD )
    call S % Deallocate_RK_C ( )
    end associate !-- ID, etc.

  end subroutine DeallocateStorage


!   subroutine Compute_1D &
!                ( S, Current_1D, Grid, Time, TimeStep, GeometryOption, &
!                  UseLimiterParameterOption, iStrgeometryValueOption )

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
!       iStrgeometryValueOption

!     integer ( KDI ) :: &
!       iF, &  !-- iField
!       iStrg     !-- iStorage
!     integer ( KDI ), dimension ( size ( Current_1D ) ) :: &
!       nValues, &
!       nEquations
!     type ( StorageForm ), dimension ( size ( Current_1D ) ) :: &
!       Solution
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     associate &
!       ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeStep ) )
!     call Timer % Start ( )

!     S % Current_1D => Current_1D
!     S % Grid => Grid

!     associate ( nStorages => size ( Current_1D ) )

!     S % Geometry => null ( )
!     if ( present ( GeometryOption ) .and. present ( iStrgeometryValueOption ) ) &
!     then
!       S % Geometry => GeometryOption
!       S % iStrgeometryValue = iStrgeometryValueOption
!     end if

!     if ( .not. allocated ( S % ApplyDivergence_1D ) ) &
!       allocate ( S % ApplyDivergence_1D ( nStorages ) )
!     if ( .not. allocated ( S % ApplySources_1D ) ) &
!       allocate ( S % ApplySources_1D ( nStorages ) )
!     if ( .not. allocated ( S % ApplyRelaxation_1D ) ) &
!       allocate ( S % ApplyRelaxation_1D ( nStorages ) )

!     !-- Allocate Solution and initialize from Current_1D

!     do iStrg = 1, nStorages
!       associate &
!         ( C   => Current_1D ( iStrg ) % Pointer )
!       associate &
!         ( iaC => C % iaConserved, &
!           nE  => C % N_CONSERVED, &  !-- nEquations
!           nV  => C % nValues )

!       call Solution ( iStrg ) % Initialize ( [ nV, nE ] )
!       do iF = 1, C % N_CONSERVED
!         associate &
!           ( CV => C % Value ( :, iaC ( iF ) ), &
!             SV => Solution ( iStrg ) % Value ( :, iF ) )
!         call Copy ( CV, SV )
!         end associate !-- CV, etc.
!       end do !-- iF

!       end associate !-- iaC, etc.
!       end associate !-- C, etc.
!     end do !-- iStrg

!     !-- Compute Solution

!     call S % SetDivergence ( UseLimiterParameterOption )
!     call S % ComputeTemplate ( Solution, Time, TimeStep )
!     call S % ClearDivergence ( )

!     !-- Copy Solution to Current_1D

!     do iStrg = 1, nStorages
!       associate ( C => Current_1D ( iStrg ) % Pointer )
!       associate ( iaC => C % iaConserved )

!       do iF = 1, C % N_CONSERVED
!         associate &
!           ( SV => Solution ( iStrg ) % Value ( :, iF ), &
!             CV => C % Value ( :, iaC ( iF ) ) )
!         call Copy ( SV, CV )
!         end associate !-- YV, etc.
!       end do !-- iF

!       if ( associated ( S % Geometry ) ) then
!         call C % ComputeFromConserved ( S % iStrgeometryValue, S % Geometry )
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
!     end do !-- iStrg

!     end associate !-- nStorages

!     nullify ( S % Grid )
!     nullify ( S % Current_1D )
!     nullify ( S % Geometry )

!     call Timer % Stop ( )
!     end associate !-- Timer

!   end subroutine Compute_1D


  subroutine Compute ( S, Time, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      Time, &
      TimeStep

    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_BF

    call Show ( 'Computing a Step_RK_C_ASC', S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerStep )
    Timer_BF => PROGRAM_HEADER % TimerPointer ( S % iTimerBoundaryFluence )

    if ( associated ( Timer ) ) call Timer % Start ( )

    if ( .not. S % Allocated ) &
      call AllocateStorage ( S )

    if ( associated ( Timer_BF ) ) call Timer_BF % Start ( )
    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      call S % ClearBoundaryFluence ( S % BoundaryFluence_CSL )
    if ( associated ( Timer_BF ) ) call Timer_BF % Stop ( )

    call S % ComputeTemplate ( Time, TimeStep )

    if ( associated ( Timer_BF ) ) call Timer_BF % Start ( )
    if ( allocated ( S % BoundaryFluence_CSL ) ) &
      call S % Current_ASC % AccumulateBoundaryTally &
             ( S % BoundaryFluence_CSL )
    if ( associated ( Timer_BF ) ) call Timer_BF % Stop ( )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine Compute


  impure elemental subroutine FinalizeTemplate_C_ASC ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    call DeallocateStorage ( S )

!    if ( allocated ( S % IncrementDamping ) ) &
!      deallocate ( S % IncrementDamping )
    if ( allocated ( S % IncrementDivergence_C ) ) &
      deallocate ( S % IncrementDivergence_C )
    if ( allocated ( S % StorageDivergence_C ) ) &
      deallocate ( S % StorageDivergence_C )

    nullify ( S % Current )
    nullify ( S % Current_ASC )

    call S % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C_ASC


  subroutine SetUseLimiter ( S, UseLimiter )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      UseLimiter

    associate ( ID => S % IncrementDivergence_C )
    call ID % SetUseLimiter ( UseLimiter )
    end associate !-- ID

  end subroutine SetUseLimiter


  subroutine InitializeTimersStage ( S, BaseLevel )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
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
             ( 'ClearIncrement', S % iTimerClearIncrement, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ComputeConstraints', S % iTimerConstraints, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ApplyDivergence', S % iTimerDivergence, &
               Level = BaseLevel + 1 )
      call S % InitializeTimersDivergence &
             ( BaseLevel + 2 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ApplySources', S % iTimerSources, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'ApplyRelaxation', S % iTimerRelaxation, &
               Level = BaseLevel + 1 )
      call PROGRAM_HEADER % AddTimer &
             ( 'GhostIncrement', S % iTimerGhost, &
               Level = BaseLevel + 1 )

  end subroutine InitializeTimersStage


  subroutine InitializeTimersDivergence ( S, BaseLevel )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( .not. allocated ( S % StorageDivergence_C ) ) &
      return

    associate ( SD => S % StorageDivergence_C )
    call SD % InitializeTimers ( BaseLevel )
    end associate !-- SD, etc.

  end subroutine InitializeTimersDivergence


  subroutine LoadSolution ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    if ( S % Current % AllocatedDevice ) &
      call S % Current % UpdateDevice ( )

    call S % LoadSolution_C ( S % Solution, S % Current )

  end subroutine LoadSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    call S % StoreSolution_C ( S % Current, S % Solution )
    
    if ( S % Current % AllocatedDevice ) &
      call S % Current % UpdateHost ( )

  end subroutine StoreSolution


  subroutine InitializeIntermediate ( S, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iStage

    if ( iStage > 1 ) &
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

    call S % ComputeStage_C_ASC ( TimeStep, iStage )

  end subroutine ComputeStage


  subroutine ComputeStage_C_ASC ( S, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    type ( TimerForm ), pointer :: &
      TimerStore, &
      TimerClear

    associate &
      ( C     => S % Current, &
        Chart => S % Chart, &
        K     => S % K ( iStage ), &
        Y     => S % Y )

    if ( iStage > 1 ) then
      TimerStore => PROGRAM_HEADER % TimerPointer &
                      ( S % iTimerStoreIntermediate )
      if ( associated ( TimerStore ) ) call TimerStore % Start ( )
      call S % StoreSolution_C ( C, Y )
      if ( associated ( TimerStore ) ) call TimerStore % Stop ( )
    end if

    TimerClear => PROGRAM_HEADER % TimerPointer ( S % iTimerClearIncrement )
    if ( associated ( TimerClear ) ) call TimerClear % Start ( )
    if ( K % AllocatedDevice ) then
      call Clear ( K % D_Selected, K % Value )
    else
      call Clear ( K % Value )
    end if
    if ( associated ( TimerClear ) ) call TimerClear % Stop ( )    

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


  subroutine LoadSolution_C ( Solution, S, Current )

    type ( StorageForm ), intent ( inout ) :: &
      Solution
    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    integer ( KDI ) :: &
      iF  !-- iField
      
    associate ( iaC => Current % iaConserved )
    do iF = 1, Current % N_CONSERVED
      
      associate &
        ( CV  => Current % Value ( :, iaC ( iF ) ), &
          D_C => Current % D_Selected ( iaC ( iF ) ), &
          SV  => Solution % Value ( :, iF ), &
          D_S => Solution % D_Selected ( iF ) )
      
      if ( Current % AllocatedDevice ) then
        call Copy ( CV, D_C, D_S, SV )
      else
        call Copy ( CV, SV )
      end if
      
      end associate !-- CV, etc.
      
    end do !-- iF
    end associate !-- iaC

  end subroutine LoadSolution_C


  subroutine StoreSolution_C ( Current, S, Solution )

    class ( CurrentTemplate ), intent ( inout ) :: &
      Current
    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( StorageForm ), intent ( in ) :: &
      Solution

    integer ( KDI ) :: &
      iF  !-- iField
    class ( GeometryFlatForm ), pointer :: &
      G

    associate ( iaC => Current % iaConserved )
    do iF = 1, Current % N_CONSERVED
      
      associate &
        ( SV  => Solution % Value ( :, iF ), &
          D_S => Solution % D_Selected ( iF ), &
          CV  => Current % Value ( :, iaC ( iF ) ), &
          D_C => Current % D_Selected ( iaC ( iF ) ) )
      
      if ( Current % AllocatedDevice ) then
        call Copy ( SV, D_S, D_C, CV )
      else
        call Copy ( SV, CV )
      end if 
      
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

    if ( associated ( S % ComputeConstraints % Pointer ) ) &
      call S % ComputeConstraints % Pointer ( S )

    nullify ( G )
    
  end subroutine StoreSolution_C


!   subroutine ComputeIncrement ( S, K, Y, Time, TimeStep, iStage, iStorage )

!     class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
!       S
!     type ( StorageForm ), dimension ( :, : ), intent ( inout ) :: &
!       K
!     type ( StorageForm ), dimension ( : ), intent ( in ) :: &
!       Y
!     real ( KDR ), intent ( in ) :: &
!       Time, &
!       TimeStep
!     integer ( KDI ), intent ( in ) :: &
!       iStage, &
!       iStorage

!     integer ( KDI ) :: &
!       iF  !-- iField
!     type ( StorageForm ), allocatable :: &
!       DC  !-- DampingCoefficient
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     associate &
!       ( Timer => PROGRAM_HEADER % Timer ( S % iTimerComputeIncrement ) )
!     call Timer % Start ( )

!     associate &
!       ( C   => S % Current_1D ( iStorage ) % Pointer, &
!         iaC => S % Current_1D ( iStorage ) % Pointer % iaConserved )

!     if ( iStage > 1 ) then

!       do iF = 1, C % N_CONSERVED
!         associate &
!           ( YV => Y ( iStorage ) % Value ( :, iF ), &
!             CV => C % Value ( :, iaC ( iF ) ) )
!         call Copy ( YV, CV )
!         end associate !-- YV, etc.
!       end do !-- iF

!       if ( associated ( S % Geometry ) ) then
!         call C % ComputeFromConserved ( S % iStrgeometryValue, S % Geometry )
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

!     associate ( KG => K ( iStorage, iStage ) )

!     call Clear ( KG % Value )

!     !-- Divergence
!     if ( associated ( S % ApplyDivergence_1D ( iStorage ) % Pointer ) ) &
!       call S % ApplyDivergence_1D ( iStorage ) % Pointer &
!              ( S, KG, C, TimeStep, iStage, iStorage )

!     !-- Other explicit sources
!     if ( associated ( S % ApplySources_1D ( iStorage ) % Pointer ) ) &
!       call S % ApplySources_1D ( iStorage ) % Pointer ( S, KG, C, TimeStep )

!     !-- Relaxation
!     if ( associated ( S % ApplyRelaxation_1D ( iStorage ) % Pointer ) ) then
!       associate ( ID => S % IncrementDamping )
!       allocate ( DC )
!       call DC % Initialize ( shape ( KG % Value ), ClearOption = .true. )
!       call S % ApplyRelaxation_1D ( iStorage ) % Pointer &
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
!       iStrg     !-- iStorage
!     integer ( KDI ), dimension ( 3 ) :: &
!       nSurface

!     associate ( nStorages => size ( S % Current_1D ) )

!     allocate ( S % UseLimiterParameter ( nStorages ) )
!     S % UseLimiterParameter = .true.
!     if ( present ( UseLimiterParameterOption ) ) &
!       S % UseLimiterParameter = UseLimiterParameterOption

!     select type ( Grid => S % Grid )
!     class is ( Chart_SLD_Form )

!       if ( allocated ( S % BoundaryFluence_CSL ) ) &
!         deallocate ( S % BoundaryFluence_CSL )
!       allocate ( S % BoundaryFluence_CSL ( nStorages ) )

!       do iStrg = 1, nStorages
!         associate &
!           ( C => Grid % Atlas % Connectivity, &
!             nDimensions => Grid % nDimensions, &
!             nConserved  => S % Current_1D ( iStrg ) % Pointer % N_CONSERVED, &
!             BF => S % BoundaryFluence_CSL ( iStrg ) )
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

!       end do !-- iStrg

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

!     end associate !-- nStorages

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
      ( SV  => S % Solution % Value, &
        D_S => S % Solution % D_Selected, &
        YV  => S % Y % Value, &
        D_Y => S % Y % D_Selected )

    if ( S % Solution % AllocatedDevice ) then
      call Copy ( SV, D_S, D_Y, YV )
    else
      call Copy ( SV, YV )
    end if

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
      ( YV  => S % Y % Value, &
        D_Y => S % Y % D_Selected, &
        KV  => S % K ( iK ) % Value, &
        D_K => S % K ( iK ) % D_Selected )
    
    if ( S % Y % AllocatedDevice ) then
      call MultiplyAdd ( YV, D_Y, D_K, KV, A )
    else
      call MultiplyAdd ( YV, KV, A )
    end if
    
    end associate !-- YV, etc.

  end subroutine IncrementIntermediate_C


  subroutine ComputeStage_C &
               ( S, ID, C, Chart, K, TimeStep, iStage, GeometryOption, & 
                 GhostExchangeOption, iStrgeometryValueOption )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID    
    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( ChartTemplate ), intent ( inout ) :: &
      Chart
    type ( StorageForm ), intent ( inout ) :: &
      K
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage
    class ( GeometryFlatForm ), intent ( in ), optional :: &
      GeometryOption
    logical ( KDL ), intent ( in ), optional :: &
      GhostExchangeOption
    integer ( KDI ), intent ( in ), optional :: &
      iStrgeometryValueOption

    logical ( KDL ) :: &
      GhostExchange
    type ( TimerForm ), pointer :: &
      TimerDivergence, &
      TimerSources, &
      TimerRelaxation, &
      TimerGhost

    !-- Divergence
    if ( associated ( S % ApplyDivergence_C ) ) then
      TimerDivergence => PROGRAM_HEADER % TimerPointer ( S % iTimerDivergence )
      if ( associated ( TimerDivergence ) ) call TimerDivergence % Start ( )
      call S % ApplyDivergence_C ( ID, K, TimeStep, iStage )
      if ( associated ( TimerDivergence ) ) call TimerDivergence % Stop ( )
    end if

    !-- Other explicit sources
    if ( associated ( S % ApplySources_C ) ) then
      TimerSources => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
      if ( associated ( TimerSources ) ) call TimerSources % Start ( )
      call S % ApplySources_C ( C % Sources, K, C, TimeStep, iStage )
      if ( associated ( TimerSources ) ) call TimerSources % Stop ( )
    end if

    !-- Relaxation
    if ( associated ( S % ApplyRelaxation_C ) ) then
      TimerRelaxation => PROGRAM_HEADER % TimerPointer ( S % iTimerRelaxation )
      if ( associated ( TimerRelaxation ) ) call TimerRelaxation % Start ( )
      call S % ApplyRelaxation_C &
             ( C, C % Sources, K, Chart, TimeStep, iStage, GeometryOption, &
               iStrgeometryValueOption )
      if ( associated ( TimerRelaxation ) ) call TimerRelaxation % Stop ( )
    end if

    if ( K % AllocatedDevice ) &
      call K % UpdateHost ( )
    
    if ( associated ( S % CoarsenSingularities ) ) &
      call S % CoarsenSingularities ( K )

    GhostExchange = .true.
    if ( present ( GhostExchangeOption ) ) &
      GhostExchange = GhostExchangeOption

    if ( GhostExchange .and. associated ( S % ApplyDivergence_C ) ) then
      select type ( Chart )
      class is ( Chart_SLD_Form )
        TimerGhost => PROGRAM_HEADER % TimerPointer ( S % iTimerGhost )
        if ( associated ( TimerGhost ) ) call TimerGhost % Start ( )
        call Chart % ExchangeGhostData ( K )
        if ( associated ( TimerGhost ) ) call TimerGhost % Stop ( )
      end select !-- Grid
    end if !-- ApplyDivergence_C
    
    if ( K % AllocatedDevice ) &
      call K % UpdateDevice ( )

  end subroutine ComputeStage_C


  subroutine IncrementSolution_C ( S, B, iS )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      B
    integer ( KDI ), intent ( in ) :: &
      iS

    associate &
      ( SV  => S % Solution % Value, &
        D_S => S % Solution % D_Selected, &
        KV  => S % K ( iS ) % Value, &
        D_K => S % K ( iS ) % D_Selected )
    
    if ( S % Solution % AllocatedDevice ) then
      call MultiplyAdd ( SV, D_S, D_K, KV, B )
    else
      call MultiplyAdd ( SV, KV, B )
    end if 
    
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
    if ( S % Current % AllocatedDevice ) &
      call S % Solution % AllocateDevice ( ) 
    
    call S % Y % Initialize ( [ nValues, nEquations ] )
    if ( S % Current % AllocatedDevice ) &
      call S % Y % AllocateDevice ( ) 
    
    do iS = 1, S % nStages
      call S % K ( iS ) % Initialize ( [ nValues, nEquations ] )
      if ( S % Current % AllocatedDevice ) &
        call S % K ( iS ) % AllocateDevice ( )
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


  subroutine AllocateStorageDivergence ( SD, C, G )

    type ( StorageDivergenceForm ), intent ( inout ) :: &
      SD
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G

    call SD % Allocate &
           ( C % AllocatedDevice, C % nVariables, C % N_CONSERVED, &
             C % N_RECONSTRUCTED, C % N_SOLVER_SPEEDS, G % nVariables, &
             C % nValues )
    if ( trim ( C % RiemannSolverType ) == 'HLLC' ) &
      call SD % Allocate_HLLC &
             ( C % AllocatedDevice, C % nVariables, C % nValues )

  end subroutine AllocateStorageDivergence


  subroutine DeallocateStorageDivergence ( SD )

    type ( StorageDivergenceForm ), intent ( inout ) :: &
      SD

    call SD % Deallocate ( )

  end subroutine DeallocateStorageDivergence


  subroutine AllocateBoundaryFluence_CSL &
               ( IncrementDivergence, CSL, nEquations, BF )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      IncrementDivergence
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    integer ( KDI ), intent ( in ) :: &
      nEquations
    type ( Real_3D_Form ), dimension ( :, : ), intent ( out ), &
      allocatable :: &
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

    call IncrementDivergence % SetBoundaryFluence ( BF )

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


  subroutine DeallocateBoundaryFluence_CSL ( ID, BF )

    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ), &
      allocatable :: &
        BF
 
    call ID % ClearBoundaryFluence ( )

    if ( allocated ( BF ) ) &
      deallocate ( BF )
    
  end subroutine DeallocateBoundaryFluence_CSL

  
  subroutine AllocateMetricDerivatives ( S, ID, nValues )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID
    integer ( KDI ), intent ( in ) :: &
      nValues

    if ( ( trim ( S % Chart % CoordinateSystem ) == 'SPHERICAL' &
           .or. trim ( S % Chart % CoordinateSystem ) == 'CYLINDRICAL' ) ) &
    then
      allocate ( S % dLogVolumeJacobian_dX ( 2 ) )
      call S % dLogVolumeJacobian_dX ( 1 ) % Initialize ( nValues )
      call S % dLogVolumeJacobian_dX ( 2 ) % Initialize ( nValues )
      call ID % SetMetricDerivatives ( S % dLogVolumeJacobian_dX )
    end if

  end subroutine AllocateMetricDerivatives


  subroutine DeallocateMetricDerivatives ( S, ID )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID

    call ID % ClearMetricDerivatives ( )

    if ( allocated ( S % dLogVolumeJacobian_dX ) ) &
      deallocate ( S % dLogVolumeJacobian_dX )

  end subroutine DeallocateMetricDerivatives


  subroutine AllocateStorage ( S )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S

    class ( GeometryFlatForm ), pointer :: &
      G

    S % Allocated = .true.

    call S % Allocate_RK_C ( )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

      G => Chart % Geometry ( )
      associate &
        ( ID => S % IncrementDivergence_C, &
          SD => S % StorageDivergence_C, &
          C  => S % Current )

      call S % AllocateStorageDivergence ( SD, C, G )

      call S % AllocateBoundaryFluence &
             ( ID, Chart, C % N_CONSERVED, S % BoundaryFluence_CSL )

      call S % AllocateMetricDerivatives ( ID, C % nValues )

      end associate !-- ID, etc.
      nullify ( G )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Step_RK_C_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateStorage', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

  end subroutine AllocateStorage


  subroutine ApplyDivergence_C ( S, ID, Increment, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    type ( IncrementDivergence_FV_Form ), intent ( inout ) :: &
      ID
    type ( StorageForm ), intent ( inout ) :: &
      Increment
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iC  !-- iConserved

    call ID % Compute ( Increment, TimeStep, Weight_RK = S % B ( iStage ) )

    !-- ID % Current is not necessarily associated until after ID % Compute
    associate ( C => ID % Current )

    if ( iStage == 1 ) &
      call Clear ( C % Sources % Value ( :, : C % N_CONSERVED ) )

    do iC = 1, C % N_CONSERVED
      call RecordDivergence &
             ( C % Sources % Value ( :, iC ), Increment % Value ( :, iC ), &
               TimeStep, Weight_RK = S % B ( iStage ) )
    end do !-- iC
!call Show ( C % Sources % Value ( :, 2 ), '>>> Divergence source' )

    end associate !-- C

  end subroutine ApplyDivergence_C


  subroutine RecordDivergence ( SDV, IDV, TimeStep, Weight_RK )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      SDV 
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      IDV 
    real ( KDR ), intent ( in ) :: &
      TimeStep, &
      Weight_RK

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( SDV )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      SDV ( iV )  =  SDV ( iV )  +  Weight_RK  *  IDV ( iV )  /  TimeStep
    end do !-- iV
    !$OMP end parallel do

  end subroutine RecordDivergence


end module Step_RK_C_ASC__Template
