!-- IntegratorHeader handles metadata of an Integrator.

module IntegratorHeader_Form

  use Basics
  use UniverseHeader_Form

  implicit none
  private

  type, public :: IntegratorHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerEvolve = 0, &
      iTimerCycle = 0, &
      iTimerNewTime = 0, &
      iTimerCheckpoint = 0, &
      iTimerTally = 0, &
      iTimerWrite = 0, &
      iCycle, &
      iCheckpoint, &
      nRampCycles, &
      nWrite, &
      nTimeStepCandidates, &
      CheckpointDisplayInterval
    real ( KDR ) :: &
      StartTime, &
      FinishTime, &
      WriteTimeInterval, &
      WriteTime, &
      Time
    real ( KDR ), dimension ( : ), allocatable :: &
      TimeStepCandidate
    logical ( KDL ) :: &
      IsCheckpointTime, &
      NoWrite, &
      WriteTimeExact
    character ( LDL ), dimension ( : ), allocatable :: &
      TimeStepLabel
    type ( MeasuredValueForm ) :: &
      TimeUnit
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    class ( UniverseHeaderForm ), pointer :: &
      Universe => null ( )
  contains
    procedure, public, pass :: &
      InitializeHeader
    final :: &
      Finalize
  end type IntegratorHeaderForm

contains


  subroutine InitializeHeader &
               ( I, U, Name, TimeUnitOption, FinishTimeOption, nWriteOption )

    class ( IntegratorHeaderForm ), intent ( inout ) :: &
      I
    class ( UniverseHeaderForm ), intent ( in ), target :: &
      U
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

    I % iCycle = 0
    I % iCheckpoint = 0
    I % nRampCycles = 100
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

    I % CheckpointDisplayInterval = 100
    call PROGRAM_HEADER % GetParameter &
           ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval' )

    I % NoWrite = .false.
    call PROGRAM_HEADER % GetParameter ( I % NoWrite, 'NoWrite' )

    I % WriteTimeExact = .false.
    call PROGRAM_HEADER % GetParameter &
           ( I % WriteTimeExact, 'WriteTimeExact' )

    if ( .not. allocated ( I % TimeStepLabel ) ) then
      allocate ( I % TimeStepLabel ( 1 ) )
      I % TimeStepLabel ( 1 ) = 'Physical'
    end if
    I % nTimeStepCandidates = size ( I % TimeStepLabel )
    allocate ( I % TimeStepCandidate ( I % nTimeStepCandidates ) )

    call Show ( I % StartTime, I % TimeUnit, 'StartTime', I % IGNORABILITY )
    call Show ( I % FinishTime, I % TimeUnit, 'FinishTime', I % IGNORABILITY )
    call Show ( I % nRampCycles, 'nRampCycles', I % IGNORABILITY )
    call Show ( I % nWrite, 'nWrite', I % IGNORABILITY )
    call Show ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval', &
                I % IGNORABILITY )

    I % Universe  =>  U

  end subroutine InitializeHeader

  
  subroutine Finalize ( I )

    type ( IntegratorHeaderForm ), intent ( inout ) :: &
      I

    if ( I % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    if ( allocated ( I % GridImageStream ) ) &
      deallocate ( I % GridImageStream )

    nullify ( I % Communicator )

  end subroutine Finalize

    
end module IntegratorHeader_Form
