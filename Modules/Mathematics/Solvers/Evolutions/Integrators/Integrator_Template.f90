!-- Integrator is a template for a time-evolved system.

module Integrator_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: IntegratorTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerAdministerCheckpoint = 0, &
      iTimerComputeCycle = 0, &
      iTimerComputeTally = 0, &
      iTimerWrite = 0, &
      iTimerComputeNewTime = 0, &
      iCycle, &
      nRampCycles, &
      nWrite
    real ( KDR ) :: &
      StartTime, &
      FinishTime, &
      WriteTimeInterval, &
      WriteTime, &
      Time
    logical ( KDL ) :: &
      IsCheckpointTime, &
      NoWrite
    type ( MeasuredValueForm ) :: &
      TimeUnit
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    class ( AtlasHeaderForm ), allocatable :: &
      PositionSpace
    class ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry_ASC
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
    procedure, private, pass :: &  !-- 2
      OpenStreams
    procedure, private, pass :: &  !-- 2
      AdministerCheckpoint
!-- FIXME: Intel compiler fails to recognize concrete overriding in
!          descendants of Integrator_C_Template. Using an empty routine.
!    procedure ( CC ), private, pass, deferred :: &  !-- 2
!      ComputeCycle
    procedure, private, pass :: &  !-- 2
      ComputeCycle
!-- See FIXME above
!    procedure ( CT ), private, pass, deferred :: &  !-- 3
!      ComputeTally
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &  !-- 3
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
!    subroutine CT ( I, ComputeChangeOption )
!      use Basics
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( inout ) :: &
!        I
!      logical ( KDL ), intent ( in ), optional :: &
!        ComputeChangeOption      
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
!    subroutine CTSL ( I, TimeStep )
!      use Basics
!      import IntegratorTemplate
!      class ( IntegratorTemplate ), intent ( in ) :: &
!        I
!      real ( KDR ), intent ( inout ) :: &
!        TimeStep
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

    I % IGNORABILITY = CONSOLE % INFO_1

    if ( I % Type == '' ) &
      I % Type = 'an Integrator' 

    I % Name = Name

    call Show ( 'Initializing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    if ( .not. allocated ( I % PositionSpace ) ) then
      call Show ( 'PositionSpace must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    else
      I % Communicator => I % PositionSpace % Communicator
    end if

    if ( .not. allocated ( I % Geometry_ASC ) ) then
      call Show ( 'Geometry must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call I % OpenStreams ( )

    I % iCycle = 0
    I % nRampCycles = 1
    call PROGRAM_HEADER % GetParameter ( I % nRampCycles, 'nRampCycles' )

    I % StartTime  = 0.0_KDR
    I % FinishTime = 1.0_KDR
    I % TimeUnit   = UNIT % IDENTITY
    if ( present ( FinishTimeOption ) ) &
      I % FinishTime = FinishTimeOption
    if ( present ( TimeUnitOption ) ) &
      I % TimeUnit = TimeUnitOption
    call PROGRAM_HEADER % GetParameter &
           ( I % StartTime, 'StartTime', InputUnitOption = I % TimeUnit )    
    call PROGRAM_HEADER % GetParameter &
           ( I % FinishTime, 'FinishTime', InputUnitOption = I % TimeUnit )

    I % nWrite = 100
    if ( present ( nWriteOption ) ) &
      I % nWrite = nWriteOption
    call PROGRAM_HEADER % GetParameter ( I % nWrite, 'nWrite' )

    I % NoWrite = .false.
    call PROGRAM_HEADER % GetParameter ( I % NoWrite, 'NoWrite' )

    call PROGRAM_HEADER % AddTimer &
           ( 'AdministerCheckpoint', I % iTimerAdministerCheckpoint )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeCycle', I % iTimerComputeCycle )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeTally', I % iTimerComputeTally )
    call PROGRAM_HEADER % AddTimer &
           ( 'Write', I % iTimerWrite )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeNewTime', I % iTimerComputeNewTime )

  end subroutine InitializeTemplate


  subroutine Evolve ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    call Show ( 'Starting evolution', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    I % Time = I % StartTime
    call I % AdministerCheckpoint ( ComputeChangeOption = .false. )

    do while ( I % Time < I % FinishTime )

      call Show ( 'Computing a cycle', I % IGNORABILITY )
      call Show ( I % Name, 'Name', I % IGNORABILITY )

      call I % ComputeCycle ( )

      call Show ( 'Cycle computed', I % IGNORABILITY )
      call Show ( I % iCycle, 'iCycle', I % IGNORABILITY )
      call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY )

      if ( I % IsCheckpointTime ) &
        call I % AdministerCheckpoint ( )

    end do !-- Time < FinishTime 

  end subroutine Evolve


  impure elemental subroutine FinalizeTemplate ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    if ( allocated ( I % MomentumSpace ) ) &
      deallocate ( I % MomentumSpace ) 
    if ( allocated ( I % Geometry_ASC ) ) &
      deallocate ( I % Geometry_ASC )
    if ( allocated ( I % PositionSpace ) ) &
      deallocate ( I % PositionSpace ) 
    if ( allocated ( I % GridImageStream ) ) &
      deallocate ( I % GridImageStream )

    nullify ( I % Communicator )

    if ( I % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

  end subroutine FinalizeTemplate


  subroutine OpenStreams ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    logical ( KDL ) :: &
      VerboseStream
    character ( LDF ) :: &
      OutputDirectory

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )

    allocate ( I % GridImageStream )
    associate ( GIS => I % GridImageStream )
    call GIS % Initialize &
           ( I % Name, CommunicatorOption = I % Communicator, &
             WorkingDirectoryOption = OutputDirectory )

    VerboseStream = .false.
    call PROGRAM_HEADER % GetParameter ( VerboseStream, 'VerboseStream' )

    associate ( iS => 1 )  !-- iStream

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call PS % OpenStream &
             ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
    class default
      call Show ( 'Atlas type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
      call Show ( 'OpenStreams', 'subroutine', CONSOLE % ERROR ) 
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

  end subroutine OpenStreams


  subroutine AdministerCheckpoint ( I, ComputeChangeOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

    real ( KDR ), dimension ( PROGRAM_HEADER % nTimers ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerAdministerCheckpoint ) )
    call Timer % Start ( )

    call I % ComputeTally ( ComputeChangeOption = ComputeChangeOption )
    if ( associated ( I % SetReference ) ) &
      call I % SetReference ( )
    call I % Write ( )
    call PROGRAM_HEADER % ShowStatistics &
           ( I % IGNORABILITY, &
             CommunicatorOption = PROGRAM_HEADER % Communicator, &
             MaxTimeOption = MaxTime, MinTimeOption = MinTime, &
             MeanTimeOption = MeanTime )
    call I % RecordTimeSeries ( MaxTime, MinTime, MeanTime )

    I % IsCheckpointTime = .false.
    if ( I % Time < I % FinishTime ) then
      call I % SetWriteTimeInterval ( )
      I % WriteTime &
        = min ( I % Time + I % WriteTimeInterval, I % FinishTime )
      call Show ( I % WriteTimeInterval, I % TimeUnit, 'WriteTimeInterval', &
                  I % IGNORABILITY )
      call Show ( I % WriteTime, I % TimeUnit, 'Next WriteTime', &
                  I % IGNORABILITY )
    else 
      call Show ( 'FinishTime reached', I % IGNORABILITY )
      call Show ( I % iCycle, 'iCycle', I % IGNORABILITY )
      call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY )
      call I % WriteTimeSeries ( )
    end if

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine AdministerCheckpoint


!-- See FIXME above
  subroutine ComputeCycle ( I )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
  end subroutine ComputeCycle


!-- See FIXME above
  subroutine ComputeTally ( I, ComputeChangeOption )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption      
  end subroutine ComputeTally


  subroutine Write ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    if ( I % NoWrite ) return

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerWrite ) )
    call Timer % Start ( )

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
        call MS % Write ( iStream = iS )
      class default
        call Show ( 'Bundle type not found', CONSOLE % ERROR )
        call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
        call Show ( 'Write', 'subroutine', CONSOLE % ERROR ) 
        call PROGRAM_HEADER % Abort ( )
      end select !-- MS
    end if !-- MomentumSpace

    end associate !-- GIS, etc.

    call Timer % Stop ( )
    end associate !-- Timer

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


  subroutine ComputeNewTime ( I, TimeNew, HoldCheckpointSolveOption )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( out ) :: &
      TimeNew
    logical ( KDL ), dimension ( : ), intent ( inout ), optional :: &
      HoldCheckpointSolveOption

    real ( KDR ) :: &
      TimeStep

    associate &
      ( Timer => PROGRAM_HEADER % Timer ( I % iTimerComputeNewTime ) )
    call Timer % Start ( )

    call Show ( 'Computing TimeNew', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

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

    if ( I % Time + TimeStep > I % WriteTime ) then
      call Show ( 'WriteTime encountered', I % IGNORABILITY ) 
!      if ( present ( HoldCheckpointSolveOption ) ) &
!        HoldCheckpointSolveOption ( 2 : C % nLevels ) = .true.
      TimeStep = I % WriteTime - I % Time
      call Show ( TimeStep, I % TimeUnit, 'Modified TimeStep', &
                  I % IGNORABILITY )
    end if

    TimeNew = I % Time + TimeStep
    call Show ( TimeNew, I % TimeUnit, 'TimeNew', I % IGNORABILITY )

!    end associate !-- C

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeNewTime


  subroutine ComputeTimeStep ( I, TimeStep )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( out ) :: &
      TimeStep

    real ( KDR ) :: &
      RampFactor
    type ( CollectiveOperation_R_Form ) :: &
      CO

    TimeStep = huge ( 0.0_KDR )

    call I % ComputeTimeStepLocal ( TimeStep )

    call CO % Initialize &
           ( I % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    CO % Outgoing % Value ( 1 ) = TimeStep

    call CO % Reduce ( REDUCTION % MIN )

    TimeStep = CO % Incoming % Value ( 1 )
    call Show ( TimeStep, I % TimeUnit, 'Physical TimeStep', I % IGNORABILITY )

    RampFactor &
      = min ( real ( I % iCycle + 1, KDR ) / I % nRampCycles, 1.0_KDR )
    if ( RampFactor < 1.0_KDR ) then
      TimeStep &
        = RampFactor * TimeStep
      call Show ( TimeStep, I % TimeUnit, 'Ramped TimeStep', I % IGNORABILITY )
    end if

  end subroutine ComputeTimeStep


!-- See FIXME above
  subroutine ComputeTimeStepLocal ( I, TimeStep )
    class ( IntegratorTemplate ), intent ( in ) :: &
      I
    real ( KDR ), intent ( inout ) :: &
      TimeStep
  end subroutine ComputeTimeStepLocal


end module Integrator_Template
