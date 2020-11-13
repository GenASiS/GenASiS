!-- TimeSeries implements a time series of memory and timing information.

module TimeSeries_Form

  use Basics
  use IntegratorHeader_Form

  implicit none
  private

  type, public :: TimeSeriesForm
    integer ( KDI ) :: &
      IGNORABILITY, &
      N_SERIES_BASIC     = 0, &
      N_SERIES_TIMER     = 0, &
      N_SERIES_TIME_STEP = 0
    integer ( KDI ) :: &
      TIME, &
      CYCLE, &
      MEMORY_MAX_HWM, &
      MEMORY_MIN_HWM, &
      MEMORY_MEAN_HWM, &
      MEMORY_MAX_RSS, &
      MEMORY_MIN_RSS, &
      MEMORY_MEAN_RSS
    integer ( KDI ) :: &
      iTime = 0
    type ( StorageForm ), allocatable :: &
      SeriesBasic, &
      SeriesTimerMax, &
      SeriesTimerMin, &
      SeriesTimerMean, &
      SeriesTimeStep
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    class ( IntegratorHeaderForm ), pointer :: &
      Integrator => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    generic, public :: &
      Initialize => InitializeBasic
    procedure, public, pass :: &
      Record
    procedure, public, pass :: &
      Write
    final :: &
      Finalize
  end type TimeSeriesForm

contains


  subroutine InitializeBasic ( TS, I )

    class ( TimeSeriesForm ), intent ( inout ) :: &
      TS
    class ( IntegratorHeaderForm ), intent ( in ), target :: &
      I

    integer ( KDI ) :: &
      iT, &  !-- iTimer
      nTimes
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    TS % IGNORABILITY = CONSOLE % INFO_2

    if ( TS % Type == '' ) &
      TS % Type = 'a TimeSeries' 

    TS % Name = 'TimeSeries_' // trim ( I % Name ) 

    call Show ( 'Initializing ' // trim ( TS % Type ), TS % IGNORABILITY )
    call Show ( TS % Name, 'Name', TS % IGNORABILITY )

    !-- Safe margin for cases where nWrite is only an estimate
    nTimes = max ( 50 * I % nWrite, 1000 )

    !-- SeriesBasic

    TS % TIME            = 1
    TS % CYCLE           = 2
    TS % MEMORY_MAX_HWM  = 3
    TS % MEMORY_MIN_HWM  = 4
    TS % MEMORY_MEAN_HWM = 5
    TS % MEMORY_MAX_RSS  = 6
    TS % MEMORY_MIN_RSS  = 7
    TS % MEMORY_MEAN_RSS = 8
    TS % N_SERIES_BASIC  = 8

    allocate ( SeriesName ( TS % N_SERIES_BASIC ) )
    SeriesName ( TS % TIME  )           = 'Time'
    SeriesName ( TS % CYCLE )           = 'Cycle'
    SeriesName ( TS % MEMORY_MAX_HWM )  = 'Memory_Max_HWM'
    SeriesName ( TS % MEMORY_MIN_HWM )  = 'Memory_Min_HWM'
    SeriesName ( TS % MEMORY_MEAN_HWM ) = 'Memory_Mean_HWM'
    SeriesName ( TS % MEMORY_MAX_RSS )  = 'Memory_Max_RSS'
    SeriesName ( TS % MEMORY_MIN_RSS )  = 'Memory_Min_RSS'
    SeriesName ( TS % MEMORY_MEAN_RSS ) = 'Memory_Mean_RSS'

    allocate ( SeriesUnit ( TS % N_SERIES_BASIC ) )
    SeriesUnit ( TS % TIME  )           = I % TimeUnit
    SeriesUnit ( TS % CYCLE )           = UNIT % IDENTITY
    SeriesUnit ( TS % MEMORY_MAX_HWM )  = UNIT % KILOBYTE
    SeriesUnit ( TS % MEMORY_MIN_HWM )  = UNIT % KILOBYTE
    SeriesUnit ( TS % MEMORY_MEAN_HWM ) = UNIT % KILOBYTE
    SeriesUnit ( TS % MEMORY_MAX_RSS )  = UNIT % KILOBYTE
    SeriesUnit ( TS % MEMORY_MIN_RSS )  = UNIT % KILOBYTE
    SeriesUnit ( TS % MEMORY_MEAN_RSS ) = UNIT % KILOBYTE

    allocate ( TS % SeriesBasic )
    associate ( SB => TS % SeriesBasic )
    call SB % Initialize &
           ( [ nTimes, TS % N_SERIES_BASIC ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Basic', &
             ClearOption = .true. )

    deallocate ( SeriesUnit )
    deallocate ( SeriesName )

    !-- SeriesTimer

    TS % N_SERIES_TIMER = PROGRAM_HEADER % nTimers

    allocate ( SeriesName ( TS % N_SERIES_TIMER ) )
    do iT = 1, TS % N_SERIES_TIMER
      SeriesName ( iT ) = PROGRAM_HEADER % Timer ( iT ) % Name
    end do !-- iT

    allocate ( SeriesUnit ( TS % N_SERIES_TIMER ) )
    do iT = 1, TS % N_SERIES_TIMER
      SeriesUnit ( iT ) = UNIT % WALL_TIME
    end do !-- iT

    allocate ( TS % SeriesTimerMax )
    allocate ( TS % SeriesTimerMin )
    allocate ( TS % SeriesTimerMean )
    allocate ( TS % SeriesTimeStep )
    associate &
      ( ST_Max  => TS % SeriesTimerMax, &
        ST_Min  => TS % SeriesTimerMin, &
        ST_Mean => TS % SeriesTimerMean, &
        ST_Step => TS % SeriesTimeStep )
    call ST_Max % Initialize &
           ( [ nTimes, TS % N_SERIES_TIMER ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Max_WallTimePerCycle', &
             ClearOption = .true. )
    call ST_Min % Initialize &
           ( [ nTimes, TS % N_SERIES_TIMER ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Min_WallTimePerCycle', &
             ClearOption = .true. )
    call ST_Mean % Initialize &
           ( [ nTimes, TS % N_SERIES_TIMER ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Mean_WallTimePerCycle', &
             ClearOption = .true. )

    deallocate ( SeriesUnit )
    deallocate ( SeriesName )

    !-- SeriesTimeStep

    TS % N_SERIES_TIME_STEP  =  I % nTimeStepCandidates

    allocate ( SeriesName ( TS % N_SERIES_TIME_STEP ) )
    allocate ( SeriesUnit ( TS % N_SERIES_TIME_STEP ) )

    SeriesName  =  I % TimeStepLabel
    SeriesUnit  =  I % TimeUnit

    call ST_Step % Initialize &
           ( [ nTimes, TS % N_SERIES_TIME_STEP ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'TimeStep', &
             ClearOption = .true. )

    deallocate ( SeriesUnit )
    deallocate ( SeriesName )

    !-- GridImageStream and Curve

    TS % Integrator => I

    if ( I % Communicator % Rank == CONSOLE % DisplayRank ) then
      allocate ( TS % GridImageStream )
      allocate ( TS % CurveImage )
      associate &
        ( GIS_I  => I  % GridImageStream, &
          GIS_TS => TS % GridImageStream, &
          CI => TS % CurveImage )
      call TS % GridImageStream % Initialize &
             ( trim ( I % Name ) // '_TimeSeries', &
               WorkingDirectoryOption = GIS_I % WorkingDirectory )
      call CI % Initialize ( GIS_TS ) 
      call CI % AddStorage ( SB )
      call CI % AddStorage ( ST_Max )
      call CI % AddStorage ( ST_Min )
      call CI % AddStorage ( ST_Mean )
      call CI % AddStorage ( ST_Step )
      end associate !-- GIS_TS, etc.
    end if !-- output rank

    !-- Cleanup

    end associate !-- SB
    end associate !-- ST_Max, etc.

  end subroutine InitializeBasic


  subroutine Record ( TS, MaxTime, MinTime, MeanTime )

    class ( TimeSeriesForm ), intent ( inout ) :: &
      TS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iT  !-- iTimer
    type ( MeasuredValueForm )  :: &
      Memory_HWM, &
      Memory_Max_HWM, &
      Memory_Min_HWM, &
      Memory_Mean_HWM, &
      Memory_RSS, &
      Memory_Max_RSS, &
      Memory_Min_RSS, &
      Memory_Mean_RSS

    associate &
      ( I   => TS % Integrator, &
        SBV => TS % SeriesBasic % Value, &
        ST_Max   => TS % SeriesTimerMax, &
        STV_Max  => TS % SeriesTimerMax % Value, &
        STV_Min  => TS % SeriesTimerMin % Value, &
        STV_Mean => TS % SeriesTimerMean % Value, &
        STSV => TS % SeriesTimeStep % Value, &
        iV  => TS % iTime )

    iV = iV + 1

    call Show ( 'Recording TimeSeries data', TS % IGNORABILITY )
    call Show ( TS % Name, 'Name', TS % IGNORABILITY )
    call Show ( iV, 'iTime', TS % IGNORABILITY )

    if ( iV > size ( SBV, dim = 1 ) ) then
      call Show ( 'Too many time series entries', CONSOLE % ERROR )
      call Show ( 'TimeSeries_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Record', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    SBV ( iV, TS % TIME )  = I % Time
    SBV ( iV, TS % CYCLE ) = I % iCycle

    call GetMemoryUsage &
           ( Memory_HWM, Memory_RSS, TS % IGNORABILITY, &
             PROGRAM_HEADER % Communicator, &
             Max_HWM_Option  = Memory_Max_HWM, &
             Min_HWM_Option  = Memory_Min_HWM, &
             Mean_HWM_Option = Memory_Mean_HWM, &
             Max_RSS_Option  = Memory_Max_RSS, &
             Min_RSS_Option  = Memory_Min_RSS, &
             Mean_RSS_Option = Memory_Mean_RSS )

    SBV ( iV, TS % MEMORY_MAX_HWM )  = Memory_Max_HWM
    SBV ( iV, TS % MEMORY_MIN_HWM )  = Memory_Min_HWM
    SBV ( iV, TS % MEMORY_MEAN_HWM ) = Memory_Mean_HWM
    SBV ( iV, TS % MEMORY_MAX_RSS )  = Memory_Max_RSS
    SBV ( iV, TS % MEMORY_MIN_RSS )  = Memory_Min_RSS
    SBV ( iV, TS % MEMORY_MEAN_RSS ) = Memory_Mean_RSS

    if ( I % iCycle > 0 ) then

      do iT = 1, TS % N_SERIES_TIMER
        STV_Max  ( iV, iT ) = MaxTime  ( iT ) / I % iCycle
        STV_Min  ( iV, iT ) = MinTime  ( iT ) / I % iCycle
        STV_Mean ( iV, iT ) = MeanTime ( iT ) / I % iCycle
        call Show ( STV_Max ( iV, iT ), &
                    trim ( ST_Max % Variable ( iT ) ) // ' per cycle', &
                    TS % IGNORABILITY )
      end do !-- iT

      STSV ( iV, : )  =  I % TimeStepCandidate

    end if

    end associate !-- I, etc.

  end subroutine Record


  subroutine Write ( TS )

    class ( TimeSeriesForm ), intent ( inout ) :: &
      TS

    if ( .not. allocated ( TS % GridImageStream ) ) &
      return

    call Show ( 'Writing a ConservationLawTimeSeries', TS % IGNORABILITY )

    associate &
      ( GIS => TS % GridImageStream, &
        SB  => TS % SeriesBasic, &
        CI => TS % CurveImage )

    call GIS % Open ( GIS % ACCESS_CREATE, SeriesOption = .false. )
    call CI % ClearGrid ( )
    call CI % SetGridWrite  &
           ( Directory = 'TimeSeries', &
             NodeCoordinate = SB % Value ( 1 : TS % iTime, TS % TIME ), &
             nProperCells = TS % iTime, oValue = 0, &
             CoordinateUnitOption = SB % Unit ( TS % TIME ), &
             CoordinateLabelOption = 't' )
    call CI % Write ( )
    call GIS % Close ( ) 

    end associate !-- GIS, etc.

  end subroutine Write


  impure elemental subroutine Finalize ( TS )

    type ( TimeSeriesForm ), intent ( inout ) :: &
      TS

    if ( allocated ( TS % CurveImage ) ) &
      deallocate ( TS % CurveImage )
    if ( allocated ( TS % GridImageStream ) ) &
      deallocate ( TS % GridImageStream )
    if ( allocated ( TS % SeriesTimeStep ) ) &
      deallocate ( TS % SeriesTimeStep )
    if ( allocated ( TS % SeriesTimerMean ) ) &
      deallocate ( TS % SeriesTimerMean )
    if ( allocated ( TS % SeriesTimerMin ) ) &
      deallocate ( TS % SeriesTimerMin )
    if ( allocated ( TS % SeriesTimerMax ) ) &
      deallocate ( TS % SeriesTimerMax )
    if ( allocated ( TS % SeriesBasic ) ) &
      deallocate ( TS % SeriesBasic )

    nullify ( TS % Integrator )

    if ( TS % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( TS % Type ), TS % IGNORABILITY )
    call Show ( TS % Name, 'Name', TS % IGNORABILITY )

  end subroutine Finalize


end module TimeSeries_Form
