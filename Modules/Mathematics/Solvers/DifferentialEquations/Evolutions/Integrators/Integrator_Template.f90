!-- Integrator is a template for a time-evolved system.

module Integrator_Template

  use Basics
  use Manifolds
  use EvolutionBasics

  implicit none
  private

  type, public, extends ( IntegratorHeaderForm ), abstract :: &
    IntegratorTemplate
      class ( AtlasHeaderForm ), allocatable :: &
        PositionSpace
      class ( BundleHeaderForm ), allocatable :: &
        MomentumSpace
      class ( TimeSeriesForm ), allocatable :: &
        TimeSeries
      procedure ( OGIS ), pointer :: &
        OpenGridImageStreams => null ( )
      procedure ( OMS ), pointer :: &
        OpenManifoldStreams => null ( )
      procedure ( A ), pointer:: &
        Analyze => null ( )
      procedure ( W ), pointer:: &
        Write => null ( )
      procedure ( R ), pointer:: &
        Read => null ( )
      procedure ( SWTI ), pointer :: &
        SetCheckpointTimeInterval => null ( )
      procedure ( CTSL ), pointer :: &
        ComputeTimeStepLocal => null ( )
      procedure ( SR ), pointer :: &
        SetReference => null ( )
      procedure ( SI ), public, pointer :: &
        SetInitial => null ( )
      procedure ( RI ), public, pointer :: &
        ResetInitial => null ( )
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate
    procedure, public, pass :: & !-- 1
      Evolve
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate
    procedure, public, pass :: &  !-- 2
      OpenGridImageStreamsTemplate
    procedure, public, pass :: &  !-- 2
      OpenManifoldStreamsTemplate
    procedure, public, pass :: &  !-- 2
      InitializeTimers
    procedure, public, pass :: &  !-- 2
      PrepareInitial
    procedure ( PE ), private, pass, deferred :: &  !-- 2
      PrepareEvolution
    procedure, private, pass :: &  !-- 2
      AdministerCheckpoint
!-- FIXME: Intel compiler fails to recognize concrete overriding in
!          descendants of Integrator_C_Template. Using an empty routine.
!    procedure ( CC ), private, pass, deferred :: &  !-- 2
!      ComputeCycle
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      InitializeStepTimers
    procedure ( UH ), private, pass, deferred :: &  !-- 3
      UpdateHost
!-- See FIXME above
!    procedure ( CT ), private, pass, deferred :: &  !-- 3
!      ComputeTally
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, public, pass :: &  !-- 3
      AnalyzeTemplate
    procedure, public, pass :: &  !-- 3
      WriteTemplate
    procedure, public, pass :: &  !-- 3
      ReadTemplate
    procedure, private, pass :: &  !-- 3
      RecordTimeSeries
    procedure, private, pass :: &  !-- 3
      WriteTimeSeries
    procedure, public, pass :: &
      PrepareCycle
    procedure, public, pass :: &  !-- 3
      ComputeNewTime
    procedure, public, pass :: &
      ComputeTimeStep
  end type IntegratorTemplate

  abstract interface 

    subroutine OGIS ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine OGIS

    subroutine OMS ( I, VerboseStreamOption )
      use Basics
      import IntegratorTemplate  
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
      logical ( KDL ), intent ( in ), optional :: &
        VerboseStreamOption
    end subroutine OMS

    subroutine A ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine A

    subroutine W ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine W

    subroutine R ( I, ReadFrom, Time, CycleNumber )
      use Basics
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
      integer ( KDI ), intent ( in ) :: &
        ReadFrom
      type ( MeasuredValueForm ), intent ( out ) :: &
        Time
      integer ( KDI ), intent ( out ) :: &
        CycleNumber
    end subroutine R

    subroutine SWTI ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine SWTI

    subroutine CTSL ( I, TimeStepCandidate )
      use Basics
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ), target :: &
        I
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        TimeStepCandidate
    end subroutine CTSL

    subroutine SR ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine SR
    
    subroutine PE ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine PE

!-- See FIXME above
!    subroutine CC ( I )
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( inout ) :: &
!        I
!    end subroutine CC

    subroutine UH ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine UH

!-- See FIXME above
!    subroutine CT ( I, ComputeChangeOption, IgnorabilityOption )
!      use Basics
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( inout ) :: &
!        I
!      logical ( KDL ), intent ( in ), optional :: &
!        ComputeChangeOption      
!      integer ( KDI ), intent ( in ), optional :: &
!        IgnorabilityOption
!    end subroutine CT

  end interface


  interface

    subroutine SI ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine SI

    subroutine RI ( I, RestartFrom, RestartTime, CycleNumber )
      use Basics
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
      integer ( KDI ), intent ( in ) :: &
        RestartFrom
      type ( MeasuredValueForm ), intent ( out ) :: &
        RestartTime
      integer ( KDI ), intent ( out ) :: &
        CycleNumber
    end subroutine RI

  end interface


    private :: &
      ResetInitial, &
      SetCheckpointTimeInterval


contains


  subroutine InitializeTemplate &
               ( I, U, Name, TimeUnitOption, FinishTimeOption, nWriteOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    class ( UniverseHeaderForm ), intent ( in ) :: &
      U
    character ( * ), intent ( in )  :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption
    
    if ( .not. allocated ( I % TimeSeries ) ) &
      allocate ( TimeSeriesForm :: I % TimeSeries )
      !-- Initialized below

    call I % InitializeHeader &
           ( U, Name, TimeUnitOption, FinishTimeOption, nWriteOption )

    if ( .not. allocated ( I % PositionSpace ) ) then
      call Show ( 'PositionSpace must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    else
      I % Communicator => I % PositionSpace % Communicator
    end if

    if ( .not. associated ( I % OpenGridImageStreams ) ) &
      I % OpenGridImageStreams => OpenGridImageStreamsTemplate
    if ( .not. associated ( I % OpenManifoldStreams ) ) &
      I % OpenManifoldStreams => OpenManifoldStreamsTemplate
    if ( .not. associated ( I % Analyze ) ) &
      I % Analyze => AnalyzeTemplate
    if ( .not. associated ( I % Write ) ) &
      I % Write => WriteTemplate
    if ( .not. associated ( I % Read ) ) &
      I % Read => ReadTemplate
    if ( .not. associated ( I % ResetInitial ) ) &
      I % ResetInitial => ResetInitial
    if ( .not. associated ( I % SetCheckpointTimeInterval ) ) &
      I % SetCheckpointTimeInterval => SetCheckpointTimeInterval

    call I % OpenGridImageStreams ( )

    !-- if allocated above, initialize
    select type ( TS => I % TimeSeries )
    type is ( TimeSeriesForm )
      call TS % Initialize ( I )
    end select !-- TS

  end subroutine InitializeTemplate


  subroutine Evolve ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    real ( KDR ) :: &
      TimeStepRatio
    type ( TimerForm ), pointer :: &
      Timer

    call I % OpenManifoldStreams ( )
    call I % InitializeTimers ( )

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerEvolve )
    if ( associated ( Timer ) ) call Timer % Start ( )   

    call I % PrepareInitial ( )
    call I % PrepareEvolution ( )
    call I % AdministerCheckpoint ( ComputeChangeOption = .false. )

    call Show ( 'Starting evolution', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    do while ( I % Time < I % FinishTime .and. I % iCycle < I % FinishCycle )
      call Show ( 'Computing a cycle', I % IGNORABILITY + 1 )

      call I % ComputeCycle ( )

      call Show ( 'Cycle computed', I % IGNORABILITY + 1 )
      call Show ( I % iCycle, 'iCycle', I % IGNORABILITY + 1 )
      call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY + 1 )

      TimeStepRatio  &
        =  minval ( I % TimeStepCandidate ) &
             / max ( I % CheckpointTimeInterval, sqrt ( tiny ( 0.0_KDR ) ) )
      if ( TimeStepRatio  <  1.0e-6  *  I % nWrite ) then
        call I % AdministerCheckpoint ( )
        call Show ( 'TimeStepRatio too small', CONSOLE % WARNING )
        call Show ( TimeStepRatio, 'TimeStepRatio', CONSOLE % WARNING )
        exit
      end if

!call I % Write ( )
      if ( I % IsCheckpointTime ) &
        call I % AdministerCheckpoint ( )

    end do !-- Time < FinishTime 

    if ( associated ( Timer ) ) call Timer % Stop ( )   

  end subroutine Evolve


  impure elemental subroutine FinalizeTemplate ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    if ( allocated ( I % MomentumSpace ) ) &
      deallocate ( I % MomentumSpace ) 
    if ( allocated ( I % PositionSpace ) ) &
      deallocate ( I % PositionSpace ) 
    if ( allocated ( I % TimeStepLabel ) ) &
      deallocate ( I % TimeStepLabel )
    if ( allocated ( I % TimeStepCandidate ) ) &
      deallocate ( I % TimeStepCandidate )

  end subroutine FinalizeTemplate


  subroutine OpenGridImageStreamsTemplate ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    character ( LDF ) :: &
      OutputDirectory

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )

    allocate ( I % GridImageStream )
    associate ( GIS => I % GridImageStream )
    call GIS % Initialize &
           ( I % Name, CommunicatorOption = I % Communicator, &
             WorkingDirectoryOption = OutputDirectory )
    end associate !-- GIS

  end subroutine OpenGridImageStreamsTemplate


  subroutine OpenManifoldStreamsTemplate ( I, VerboseStreamOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      VerboseStreamOption

    logical ( KDL ) :: &
      VerboseStream

    associate ( GIS => I % GridImageStream )

    VerboseStream = .false.
    if ( present ( VerboseStreamOption ) ) &
      VerboseStream = VerboseStreamOption
    call PROGRAM_HEADER % GetParameter ( VerboseStream, 'VerboseStream' )

    associate ( iS => 1 )  !-- iStream

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call PS % OpenStream &
             ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
    class default
      call Show ( 'Atlas type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'OpenManifoldStreams', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- PS

    if ( allocated ( I % MomentumSpace ) ) then
      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )
        call MS % OpenStream ( GIS, iStream = iS )
      class default
        call Show ( 'Bundle type not found', CONSOLE % ERROR )
        call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
        call Show ( 'OpenStreams', 'subroutine', CONSOLE % ERROR ) 
        call PROGRAM_HEADER % Abort ( )
      end select !-- MS
    end if !-- MomentumSpace

    end associate !-- iS
    end associate !-- GIS

  end subroutine OpenManifoldStreamsTemplate


  subroutine InitializeTimers ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      BaseLevel

    BaseLevel = 0

    if ( I % iTimerEvolve > 0 ) &
      return

    call PROGRAM_HEADER % AddTimer &
           ( 'Evolve', I % iTimerEvolve, Level = BaseLevel )
      call PROGRAM_HEADER % AddTimer &
             ( 'Cycle', I % iTimerCycle, &
               Level = BaseLevel + 1 )
        call PROGRAM_HEADER % AddTimer &
             ( 'NewTime', I % iTimerNewTime, &
               Level = BaseLevel + 2 )
        call I % InitializeStepTimers ( BaseLevel + 2 )
      call PROGRAM_HEADER % AddTimer &
             ( 'Checkpoint', I % iTimerCheckpoint, &
               Level = BaseLevel + 1 )
        call PROGRAM_HEADER % AddTimer &
               ( 'Tally', I % iTimerTally, &
                 Level = BaseLevel + 2 )
        call PROGRAM_HEADER % AddTimer &
               ( 'Analyze', I % iTimerAnalyze, &
                 Level = BaseLevel + 2 )
        call PROGRAM_HEADER % AddTimer &
               ( 'Write', I % iTimerWrite, &
                 Level = BaseLevel + 2 )
        call PROGRAM_HEADER % AddTimer &
               ( 'WriteSeries', I % iTimerWriteSeries, &
                 Level = BaseLevel + 2 )

  end subroutine InitializeTimers


  subroutine PrepareInitial ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      RestartFrom, &
      CycleNumber
    type ( MeasuredValueForm ) :: &
      RestartTime

    I % Start    =  .true.
    I % Time     =  I % StartTime
    if ( associated ( I % SetInitial ) ) &
      call I % SetInitial ( )

    I % Restart  =  .false.
    RestartFrom  =  - huge ( 1 )
    call PROGRAM_HEADER % GetParameter ( RestartFrom, 'RestartFrom' )

    if ( RestartFrom >= 0 ) then
      call I % ResetInitial ( RestartFrom, RestartTime, CycleNumber )
      I % Start        =  .false.
      I % Restart      =  .true.
      I % iCheckpoint  =  RestartFrom
      I % iCycle       =  CycleNumber
      I % Time         =  RestartTime
    end if

  end subroutine PrepareInitial


  subroutine AdministerCheckpoint ( I, ComputeChangeOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

    integer ( KDI ) :: &
      iTSC, &
      TallyIgnorability, &
      StatisticsIgnorability
    real ( KDR ), dimension ( PROGRAM_HEADER % nTimers ) :: &
      MaxTime, &
      MinTime, &
      MeanTime
    logical ( KDL ) :: &
      WriteSeries
    type ( TimerForm ), pointer :: &
      Timer, &
      Timer_T, &
      Timer_A, &
      Timer_W, &
      Timer_WS
    
    Timer    => PROGRAM_HEADER % TimerPointer ( I % iTimerCheckpoint )
    Timer_T  => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    Timer_A  => PROGRAM_HEADER % TimerPointer ( I % iTimerAnalyze )
    Timer_W  => PROGRAM_HEADER % TimerPointer ( I % iTimerWrite )
    Timer_WS => PROGRAM_HEADER % TimerPointer ( I % iTimerWriteSeries )

    if ( associated ( Timer ) ) call Timer % Start ( )   

    call Show ( 'Checkpoint reached', I % IGNORABILITY )
    call Show ( I % iCheckpoint, 'iCheckpoint', I % IGNORABILITY )
    call Show ( I % iCycle, 'iCycle', I % IGNORABILITY )
    call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY )
    if ( .not. I % Start .and. .not. I % Restart ) then
      do iTSC = 1, I % nTimeStepCandidates
        call Show ( I % TimeStepCandidate ( iTSC ), I % TimeUnit, &
                    trim ( I % TimeStepLabel ( iTSC ) ) // ' TimeStep', &
                    I % IGNORABILITY )
      end do !-- iTSC
    end if

    call I % UpdateHost ( )

    if ( .not. I % Start .and. .not. I % Restart &
         .and. I % Time < I % FinishTime &
         .and. mod ( I % iCheckpoint, I % CheckpointDisplayInterval ) > 0 ) &
    then
      TallyIgnorability       =  I % IGNORABILITY + 2
      StatisticsIgnorability  =  I % IGNORABILITY + 2
      WriteSeries = .false.
    else
      TallyIgnorability       =  CONSOLE % INFO_1
      StatisticsIgnorability  =  CONSOLE % INFO_1
      WriteSeries = .true.
    end if

    if ( associated ( Timer_T ) ) call Timer_T % Start ( )   
    call I % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption  = TallyIgnorability )
    if ( associated ( Timer_T ) ) call Timer_T % Stop ( )   

    if ( associated ( Timer_A ) ) call Timer_A % Start ( )   
    call I % Analyze ( )
    if ( associated ( Timer_A ) ) call Timer_A % Stop ( )   

    if ( associated ( Timer_W ) ) call Timer_W % Start ( )   
    if ( .not. I % NoWrite ) &
      call I % Write ( )
    if ( associated ( Timer_W ) ) call Timer_W % Stop ( )   

    call PROGRAM_HEADER % ShowStatistics &
           ( StatisticsIgnorability, &
             CommunicatorOption = PROGRAM_HEADER % Communicator, &
             MaxTimeOption = MaxTime, MinTimeOption = MinTime, &
             MeanTimeOption = MeanTime )

    call I % RecordTimeSeries ( MaxTime, MinTime, MeanTime )

    if ( associated ( Timer_WS ) ) call Timer_WS % Start ( )   
    if ( WriteSeries .and. .not. I % NoWrite ) &
      call I % WriteTimeSeries ( )
    if ( associated ( Timer_WS ) ) call Timer_WS % Stop ( )   

    I % IsCheckpointTime = .false.
    if ( I % Time < I % FinishTime ) then
      call I % SetCheckpointTimeInterval ( )
      I % CheckpointTime &
        = min ( I % Time + I % CheckpointTimeInterval, I % FinishTime )
      if ( I % CheckpointTime == I % FinishTime ) &
        I % CheckpointTimeExact = .true.
      call Show ( I % CheckpointTimeInterval, I % TimeUnit, &
                  'CheckpointTimeInterval', &
                  I % IGNORABILITY )
      call Show ( I % CheckpointTime, I % TimeUnit, 'Next CheckpointTime', &
                  I % IGNORABILITY + 1 )
    else 
      call Show ( 'FinishTime reached', I % IGNORABILITY )
    end if

    I % iCheckpoint = I % iCheckpoint + 1
    
    I % Start    =  .false.
    I % Restart  =  .false.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine AdministerCheckpoint


!-- See FIXME above
  subroutine ComputeCycle ( I )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
  end subroutine ComputeCycle


  subroutine InitializeStepTimers ( I, BaseLevel )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      BaseLevel
  end subroutine InitializeStepTimers


!-- See FIXME above
  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption      
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
  end subroutine ComputeTally


  subroutine AnalyzeTemplate ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    if ( associated ( I % SetReference ) ) &
      call I % SetReference ( )

  end subroutine AnalyzeTemplate


  subroutine WriteTemplate ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    if ( allocated ( I % MomentumSpace ) ) then
      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )
        call MS % MarkFibersWritten ( )
      end select !-- MS
    end if !-- MomentumSpace

    associate &
      ( GIS => I % GridImageStream, &
        iS  => 1 )  !-- iStream
    call GIS % Open ( GIS % ACCESS_CREATE )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call PS % Write &
             ( iStream = iS, TimeOption = I % Time / I % TimeUnit, &
               CycleNumberOption = I % iCycle )
    class default
      call Show ( 'Atlas type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'Write', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- PS

    !-- Base's GIS must be closed before call to Bundle % Write ( ).
    call GIS % Close ( )

    if ( allocated ( I % MomentumSpace ) ) then
      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )
        call MS % Write &
               ( iStream = iS, TimeOption = I % Time / I % TimeUnit, &
                 CycleNumberOption = I % iCycle )
      class default
        call Show ( 'Bundle type not found', CONSOLE % ERROR )
        call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
        call Show ( 'Write', 'subroutine', CONSOLE % ERROR ) 
        call PROGRAM_HEADER % Abort ( )
      end select !-- MS
    end if !-- MomentumSpace

    end associate !-- GIS, etc.

  end subroutine WriteTemplate


  subroutine ReadTemplate ( I, ReadFrom, Time, CycleNumber )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      ReadFrom
    type ( MeasuredValueForm ), intent ( out ) :: &
      Time
    integer ( KDI ), intent ( out ) :: &
      CycleNumber

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % MarkFibersWritten ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    associate &
      ( GIS => I % GridImageStream, &
        iS  => 1 )  !-- iStream
    call GIS % Open ( GIS % ACCESS_READ, NumberOption = ReadFrom )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call PS % Read &
             ( iStream = iS, TimeOption = Time, &
               CycleNumberOption = CycleNumber )
      Time  =  Time * I % TimeUnit
    class default
      call Show ( 'Atlas type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'Read', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- PS

    !-- Base's GIS must be closed before call to Bundle % Write ( ).
    call GIS % Close ( )

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % Write &
    !            ( iStream = iS, TimeOption = I % Time / I % TimeUnit, &
    !              CycleNumberOption = I % iCycle )
    !   class default
    !     call Show ( 'Bundle type not found', CONSOLE % ERROR )
    !     call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
    !     call Show ( 'Write', 'subroutine', CONSOLE % ERROR ) 
    !     call PROGRAM_HEADER % Abort ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    end associate !-- GIS, etc.

  end subroutine ReadTemplate


  subroutine RecordTimeSeries ( I, MaxTime, MinTime, MeanTime )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iT  !-- iTimer
    real ( KDR ) :: &
      ReconstructionImbalance

    if ( .not. allocated ( I % TimeSeries ) ) &
      return

    call I % TimeSeries % Record ( MaxTime, MinTime, MeanTime )

  end subroutine RecordTimeSeries


  subroutine WriteTimeSeries ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( I % TimeSeries ) ) &
      return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerWriteSeries )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call I % TimeSeries % Write ( )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine WriteTimeSeries


  subroutine PrepareCycle ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

  end subroutine PrepareCycle


  subroutine ComputeNewTime ( I, TimeNew, HoldCheckpointSolveOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( out ) :: &
      TimeNew
    logical ( KDL ), dimension ( : ), intent ( inout ), optional :: &
      HoldCheckpointSolveOption

    real ( KDR ) :: &
      TimeStep
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerNewTime )
    if ( associated ( Timer ) ) call Timer % Start ( )

!    call Show ( 'Computing TimeNew', I % IGNORABILITY )
!    call Show ( I % Name, 'Name', I % IGNORABILITY )

    ! associate ( CFC => I % ConservedFields % Chart ( 1 ) % Element )
    ! select type ( C => I % Atlas % Chart ( 1 ) % Element )
    ! class is ( Chart_SL_Template )
    !   call I % ComputeTimeStep_CSL ( CFC, C, TimeStep )
    ! ! class is ( ChartMultiLevelTemplate )
    ! !   call I % ComputeTimeStep_CML ( CFC, C, TimeStep )
    ! end select !-- C
    ! end associate !-- CFC

    call I % ComputeTimeStep ( TimeStep )

!    associate ( C => I % Atlas % Chart ( 1 ) % Element )

    if ( I % CheckpointTimeExact ) then
      if ( I % Time + TimeStep > I % CheckpointTime ) then
        call Show ( 'CheckpointTime encountered', I % IGNORABILITY ) 
!        if ( present ( HoldCheckpointSolveOption ) ) &
!          HoldCheckpointSolveOption ( 2 : C % nLevels ) = .true.
        TimeStep = I % CheckpointTime - I % Time
        call Show ( TimeStep, I % TimeUnit, 'Modified TimeStep', &
                    I % IGNORABILITY + 1 )
      end if
    end if

    TimeNew = I % Time + TimeStep
    call Show ( TimeNew, I % TimeUnit, 'TimeNew', I % IGNORABILITY + 1 )

!    end associate !-- C

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeNewTime


  subroutine ComputeTimeStep ( I, TimeStep )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( out ) :: &
      TimeStep

    integer ( KDI ) :: &
      iTSC  !-- iTimeStepCandidate
    real ( KDR ) :: &
      RampFactor
    type ( CollectiveOperation_R_Form ) :: &
      CO

    I % TimeStepCandidate = huge ( 0.0_KDR )

    call I % ComputeTimeStepLocal ( I % TimeStepCandidate )

    call CO % Initialize &
           ( I % Communicator, nOutgoing = [ I % nTimeStepCandidates ], &
             nIncoming = [ I % nTimeStepCandidates ] )
    CO % Outgoing % Value = I % TimeStepCandidate

    call CO % Reduce ( REDUCTION % MIN )

    I % TimeStepCandidate = CO % Incoming % Value
    do iTSC = 1, I % nTimeStepCandidates
      call Show ( I % TimeStepCandidate ( iTSC ), I % TimeUnit, &
                  trim ( I % TimeStepLabel ( iTSC ) ) // ' TimeStep', &
                  I % IGNORABILITY + 1 )
    end do !-- iTSC

    TimeStep = minval ( I % TimeStepCandidate )

    RampFactor &
      = min ( real ( I % iCycle + 1, KDR ) / I % nRampCycles, 1.0_KDR )
    if ( RampFactor < 1.0_KDR ) then
      TimeStep &
        = RampFactor * TimeStep
      call Show ( TimeStep, I % TimeUnit, 'Ramped TimeStep', &
                  I % IGNORABILITY + 1 )
    end if

  end subroutine ComputeTimeStep


  subroutine ResetInitial ( I, RestartFrom, RestartTime, CycleNumber )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      RestartFrom
    type ( MeasuredValueForm ), intent ( out ) :: &
      RestartTime
    integer ( KDI ), intent ( out ) :: &
      CycleNumber

    call I % Read ( RestartFrom, RestartTime, CycleNumber )

  end subroutine ResetInitial


  subroutine SetCheckpointTimeInterval ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    I % CheckpointTimeInterval &
      = ( I % FinishTime - I % StartTime ) / I % nWrite

  end subroutine SetCheckpointTimeInterval


end module Integrator_Template
