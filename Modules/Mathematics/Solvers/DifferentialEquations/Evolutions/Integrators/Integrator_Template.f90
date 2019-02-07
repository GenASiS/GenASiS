!-- Integrator is a template for a time-evolved system.

module Integrator_Template

  use Basics
  use Manifolds
  use EvolutionBasics

  implicit none
  private

  type, public, extends ( IntegratorHeaderForm ), abstract :: &
    IntegratorTemplate
      type ( CommunicatorForm ), pointer :: &
        Communicator => null ( )
      type ( GridImageStreamForm ), allocatable :: &
        GridImageStream
      class ( AtlasHeaderForm ), allocatable :: &
        PositionSpace
      class ( BundleHeaderForm ), allocatable :: &
        MomentumSpace
      procedure ( SR ), pointer :: &
        SetReference => null ( )
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
      OpenGridImageStreams
    procedure, public, pass :: &  !-- 2
      OpenManifoldStreamsTemplate
    procedure, public, pass :: &  !-- 2
      OpenManifoldStreams
    procedure, public, pass :: &  !-- 2
      InitializeTimers
    procedure, private, pass :: &  !-- 2
      ComputeConstraints
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
!-- See FIXME above
!    procedure ( CT ), private, pass, deferred :: &  !-- 3
!      ComputeTally
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, public, pass :: &  !-- 3
      WriteTemplate
    procedure, public, pass :: &  !-- 3
      Write
!-- See FIXME above
!    procedure ( RTS ), private, pass, deferred :: &  !-- 3
!      RecordTimeSeries
    procedure, private, pass :: &  !-- 3
      RecordTimeSeries
!-- See FIXME above
!    procedure ( CC ), private, pass, deferred :: &  !-- 3
!      WriteTimeSeries
    procedure, private, pass :: &  !-- 3
      WriteTimeSeries
    procedure, private, pass :: &  !-- 3
      SetWriteTimeInterval
    procedure, public, pass :: &
      PrepareCycle
    procedure, public, pass :: &  !-- 3
      ComputeNewTime
    procedure, public, pass :: &
      ComputeTimeStep
!-- See FIXME above
!    procedure ( CTSL ), private, pass, deferred :: &
!      ComputeTimeStepLocal
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type IntegratorTemplate

  abstract interface 

    subroutine SR ( I )
      import IntegratorTemplate
      class ( IntegratorTemplate ), intent ( inout ) :: &
        I
    end subroutine SR

!-- See FIXME above
!    subroutine CC ( I )
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( inout ) :: &
!        I
!    end subroutine CC

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

!-- See FIXME above
!    subroutine RTS ( I, MaxTime, MinTime, MeanTime )
!      use Basics
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( inout ) :: &
!        I
!      real ( KDR ), dimension ( : ), intent ( in ) :: &
!        MaxTime, &
!        MinTime, &
!        MeanTime
!    end subroutine RTS

!-- See FIXME above
!    subroutine CTSL ( I, TimeStepCandidate )
!      use Basics
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( in ), target :: &
!        I
!      real ( KDR ), dimension ( : ), intent ( inout ) :: &
!        TimeStepCandidate
!    end subroutine CTSL

  end interface


contains


  subroutine InitializeTemplate &
               ( I, Name, TimeUnitOption, FinishTimeOption, nWriteOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    call I % IntegratorHeaderForm % InitializeHeader &
           ( Name, TimeUnitOption, FinishTimeOption, nWriteOption )

    if ( .not. allocated ( I % PositionSpace ) ) then
      call Show ( 'PositionSpace must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    else
      I % Communicator => I % PositionSpace % Communicator
    end if

    call I % OpenGridImageStreams ( )

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

    call Show ( 'Starting evolution', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    I % Time = I % StartTime
    call I % ComputeConstraints ( )
    call I % AdministerCheckpoint ( ComputeChangeOption = .false. )

    do while ( I % Time < I % FinishTime )
      call Show ( 'Computing a cycle', I % IGNORABILITY + 1 )

      call I % ComputeCycle ( )

      call Show ( 'Cycle computed', I % IGNORABILITY + 1 )
      call Show ( I % iCycle, 'iCycle', I % IGNORABILITY + 1 )
      call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY + 1 )

      TimeStepRatio  &
        =  minval ( I % TimeStepCandidate ) &
             / max ( I % WriteTimeInterval, sqrt ( tiny ( 0.0_KDR ) ) )
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
    if ( allocated ( I % GridImageStream ) ) &
      deallocate ( I % GridImageStream )
    if ( allocated ( I % TimeStepLabel ) ) &
      deallocate ( I % TimeStepLabel )
    if ( allocated ( I % TimeStepCandidate ) ) &
      deallocate ( I % TimeStepCandidate )

    nullify ( I % Communicator )

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


  subroutine OpenGridImageStreams ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    call I % OpenGridImageStreamsTemplate ( )

  end subroutine OpenGridImageStreams


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


  subroutine OpenManifoldStreams ( I, VerboseStreamOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      VerboseStreamOption

    call I % OpenManifoldStreamsTemplate ( VerboseStreamOption )

  end subroutine OpenManifoldStreams


  subroutine InitializeTimers ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      BaseLevel

    BaseLevel = 0

    if ( I % iTimerEvolve > 0  &
         .or.  BaseLevel > PROGRAM_HEADER % TimerLevel ) &
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
               ( 'Write', I % iTimerWrite, &
                 Level = BaseLevel + 2 )

  end subroutine InitializeTimers


  subroutine ComputeConstraints ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

  end subroutine ComputeConstraints


  subroutine AdministerCheckpoint ( I, ComputeChangeOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

    integer ( KDI ) :: &
      iTSC, &
      TallyIgnorability, &
      StatisticsIgnorability
    logical ( KDL ) :: &
      WriteTimeSeries
    real ( KDR ), dimension ( PROGRAM_HEADER % nTimers ) :: &
      MaxTime, &
      MinTime, &
      MeanTime
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerCheckpoint )
    if ( associated ( Timer ) ) call Timer % Start ( )   

    call Show ( 'Checkpoint reached', I % IGNORABILITY )
    call Show ( I % iCheckpoint, 'iCheckpoint', I % IGNORABILITY )
    call Show ( I % iCycle, 'iCycle', I % IGNORABILITY )
    call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY )
    do iTSC = 1, I % nTimeStepCandidates
      call Show ( I % TimeStepCandidate ( iTSC ), I % TimeUnit, &
                  trim ( I % TimeStepLabel ( iTSC ) ) // ' TimeStep', &
                  I % IGNORABILITY )
    end do !-- iTSC

    if ( I % Time > I % StartTime .and. I % Time < I % FinishTime &
         .and. mod ( I % iCheckpoint, I % CheckpointDisplayInterval ) > 0 ) &
    then
      TallyIgnorability      = I % IGNORABILITY + 2
      StatisticsIgnorability = I % IGNORABILITY + 2
      WriteTimeSeries        = .false.
    else
      TallyIgnorability      = CONSOLE % INFO_1
      StatisticsIgnorability = CONSOLE % INFO_1
      WriteTimeSeries        = .true.
    end if
    call I % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption  = TallyIgnorability )

    if ( associated ( I % SetReference ) ) &
      call I % SetReference ( )
    call I % Write ( )
    call PROGRAM_HEADER % ShowStatistics &
           ( StatisticsIgnorability, &
             CommunicatorOption = PROGRAM_HEADER % Communicator, &
             MaxTimeOption = MaxTime, MinTimeOption = MinTime, &
             MeanTimeOption = MeanTime )
    call I % RecordTimeSeries ( MaxTime, MinTime, MeanTime )
    if ( WriteTimeSeries ) &
      call I % WriteTimeSeries ( )

    I % IsCheckpointTime = .false.
    if ( I % Time < I % FinishTime ) then
      call I % SetWriteTimeInterval ( )
      I % WriteTime &
        = min ( I % Time + I % WriteTimeInterval, I % FinishTime )
      if ( I % WriteTime == I % FinishTime ) &
        I % WriteTimeExact = .true.
      call Show ( I % WriteTimeInterval, I % TimeUnit, 'WriteTimeInterval', &
                  I % IGNORABILITY )
      call Show ( I % WriteTime, I % TimeUnit, 'Next WriteTime', &
                  I % IGNORABILITY + 1 )
    else 
      call Show ( 'FinishTime reached', I % IGNORABILITY )
    end if

    I % iCheckpoint = I % iCheckpoint + 1

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


  subroutine WriteTemplate ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    type ( TimerForm ), pointer :: &
      Timer

    if ( I % NoWrite ) return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerWrite )
    if ( associated ( Timer ) ) call Timer % Start ( )

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

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine WriteTemplate


  subroutine Write ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    call I % WriteTemplate ( )

  end subroutine Write


!-- See FIXME above
  subroutine RecordTimeSeries ( I, MaxTime, MinTime, MeanTime )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime
  end subroutine RecordTimeSeries


!-- See FIXME above
  subroutine WriteTimeSeries ( I )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
  end subroutine WriteTimeSeries


  subroutine SetWriteTimeInterval ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    I % WriteTimeInterval &
      = ( I % FinishTime - I % StartTime ) / I % nWrite

  end subroutine SetWriteTimeInterval


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

    if ( I % WriteTimeExact ) then
      if ( I % Time + TimeStep > I % WriteTime ) then
        call Show ( 'WriteTime encountered', I % IGNORABILITY ) 
!        if ( present ( HoldCheckpointSolveOption ) ) &
!          HoldCheckpointSolveOption ( 2 : C % nLevels ) = .true.
        TimeStep = I % WriteTime - I % Time
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


!-- See FIXME above
  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )
    class ( IntegratorTemplate ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate
  end subroutine ComputeTimeStepLocal


end module Integrator_Template
