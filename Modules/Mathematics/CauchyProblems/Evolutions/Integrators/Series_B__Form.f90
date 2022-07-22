module Series_B__Form

  !-- Series_Basic__Form

  use Basics

  implicit none
  private

  type, public :: Series_B_Form
    integer ( KDI ) :: &
      IGNORABILITY, &
      N_SERIES_B     = 0, &
      N_SERIES_TIMER = 0, &
      N_SERIES_DT    = 0
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
      iTimerRecord = 0, &
      iTimerWrite  = 0, &
      iTimerRead   = 0
    integer ( KDI ) :: &
      iRecord = 0
    integer ( KDI ), pointer :: &
      iCycle => null ( )
    real ( KDR ), pointer :: &
      T => null ( )
    real ( KDR ), dimension ( : ), pointer :: &
      dT_Candidate => null ( )
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), allocatable :: &
      Basic, &
      TimerMax, &
      TimerMin, &
      TimerMean, &
      dT
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( CurveImageForm ), allocatable :: &
      CurveImage
  contains
    procedure, private, pass :: &
      Initialize_B
    generic, public :: &
      Initialize => Initialize_B
    procedure, public, pass :: &
      TimerRecord
    procedure, public, pass :: &
      TimerWrite
    procedure, public, pass :: &
      TimerRead
    procedure, public, pass :: &
      Record
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    procedure, public, pass :: &
      Restore
    final :: &
      Finalize
  end type Series_B_Form


contains


  subroutine Initialize_B &
               ( S, GIS, dT_Label, Unit_T, dT_Candidate, T, CommunicatorRank, &
                 nWrite, iCycle )

    class ( Series_B_Form ), intent ( inout ) :: &
      S
    type ( GridImageStreamForm ), intent ( in ) :: &
      GIS
    character ( * ), dimension ( : ), intent ( in ) :: &
      dT_Label
    type ( QuantityForm ), intent ( in ) :: &
      Unit_T
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      dT_Candidate
    real ( KDR ), intent ( in ), target :: &
      T
    integer ( KDI ), intent ( in ) :: &
      CommunicatorRank, &
      nWrite
    integer ( KDI ), intent ( in ), target :: &
      iCycle

    integer ( KDI ) :: &
      iT, &  !-- iTimer
      nTimes
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    S % IGNORABILITY  =  CONSOLE % INFO_1

    if ( S % Type == '' ) &
      S % Type = 'a Series_B' 

    S % Name  =  'Series'

    call Show ( 'Initializing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    !-- Safe margin for cases where nWrite is only an estimate
    nTimes  =  max ( 50 * nWrite, 1000 )

    !-- SeriesBasic

    S % TIME             =  1
    S % CYCLE            =  2
    S % MEMORY_MAX_HWM   =  3
    S % MEMORY_MIN_HWM   =  4
    S % MEMORY_MEAN_HWM  =  5
    S % MEMORY_MAX_RSS   =  6
    S % MEMORY_MIN_RSS   =  7
    S % MEMORY_MEAN_RSS  =  8
    S % N_SERIES_B       =  8

    allocate ( SeriesName ( S % N_SERIES_B ) )
    SeriesName ( S % TIME  )            =  'Time'
    SeriesName ( S % CYCLE )            =  'Cycle'
    SeriesName ( S % MEMORY_MAX_HWM )   =  'Memory_Max_HWM'
    SeriesName ( S % MEMORY_MIN_HWM )   =  'Memory_Min_HWM'
    SeriesName ( S % MEMORY_MEAN_HWM )  =  'Memory_Mean_HWM'
    SeriesName ( S % MEMORY_MAX_RSS )   =  'Memory_Max_RSS'
    SeriesName ( S % MEMORY_MIN_RSS )   =  'Memory_Min_RSS'
    SeriesName ( S % MEMORY_MEAN_RSS )  =  'Memory_Mean_RSS'

    allocate ( SeriesUnit ( S % N_SERIES_B ) )
    SeriesUnit ( S % TIME  )            =  Unit_T
    SeriesUnit ( S % CYCLE )            =  UNIT % IDENTITY
    SeriesUnit ( S % MEMORY_MAX_HWM )   =  UNIT % KILOBYTE
    SeriesUnit ( S % MEMORY_MIN_HWM )   =  UNIT % KILOBYTE
    SeriesUnit ( S % MEMORY_MEAN_HWM )  =  UNIT % KILOBYTE
    SeriesUnit ( S % MEMORY_MAX_RSS )   =  UNIT % KILOBYTE
    SeriesUnit ( S % MEMORY_MIN_RSS )   =  UNIT % KILOBYTE
    SeriesUnit ( S % MEMORY_MEAN_RSS )  =  UNIT % KILOBYTE

    allocate ( S % Basic )
    associate ( B  =>  S % Basic )
    call B % Initialize &
           ( [ nTimes, S % N_SERIES_B ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Basic', &
             ClearOption = .true. )

    deallocate ( SeriesUnit )
    deallocate ( SeriesName )

    !-- Timer

    S % N_SERIES_TIMER  =  PROGRAM_HEADER % Timer_1D % nTimers

    allocate ( SeriesName ( S % N_SERIES_TIMER ) )
    do iT  =  1, S % N_SERIES_TIMER
      SeriesName ( iT )  =  PROGRAM_HEADER % Timer_1D % Element ( iT ) % Name
    end do !-- iT

    allocate ( SeriesUnit ( S % N_SERIES_TIMER ) )
    do iT  =  1, S % N_SERIES_TIMER
      SeriesUnit ( iT )  =  UNIT % WALL_TIME
    end do !-- iT

    allocate ( S % TimerMax )
    allocate ( S % TimerMin )
    allocate ( S % TimerMean )
    associate &
      ( T_Max  => S % TimerMax, &
        T_Min  => S % TimerMin, &
        T_Mean => S % TimerMean )
    call T_Max % Initialize &
           ( [ nTimes, S % N_SERIES_TIMER ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Max_WallTimePerCycle', &
             ClearOption = .true. )
    call T_Min % Initialize &
           ( [ nTimes, S % N_SERIES_TIMER ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Min_WallTimePerCycle', &
             ClearOption = .true. )
    call T_Mean % Initialize &
           ( [ nTimes, S % N_SERIES_TIMER ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'Mean_WallTimePerCycle', &
             ClearOption = .true. )

    deallocate ( SeriesUnit )
    deallocate ( SeriesName )

    !-- SeriesTimeStep

    S % N_SERIES_DT  =  size ( dT_Label )

    allocate ( SeriesName ( S % N_SERIES_DT ) )
    allocate ( SeriesUnit ( S % N_SERIES_DT ) )

    SeriesName  =  dT_Label
    SeriesUnit  =  Unit_T

    allocate ( S % dT )
    associate ( dT => S % dT )

    call dT % Initialize &
           ( [ nTimes, S % N_SERIES_DT ], VariableOption = SeriesName, &
             UnitOption = SeriesUnit, NameOption = 'dT', &
             ClearOption = .true. )

    deallocate ( SeriesUnit )
    deallocate ( SeriesName )

    !-- GridImageStream and Curve

    if ( CommunicatorRank  ==  CONSOLE % DisplayRank ) then
      allocate ( S % GridImageStream )
      allocate ( S % CurveImage )
      associate &
        ( CI     =>  S % CurveImage, &
          GIS_S  =>  S % GridImageStream )
      call GIS_S % Initialize &
             ( trim ( GIS % Name ) // '_' // trim ( S % Name ), &
               WorkingDirectoryOption = GIS % WorkingDirectory )
      call CI % Initialize ( GIS_S ) 
      call CI % AddStorage ( B )
      call CI % AddStorage ( T_Max )
      call CI % AddStorage ( T_Min )
      call CI % AddStorage ( T_Mean )
      call CI % AddStorage ( dT )
      end associate !-- CI, etc.
    end if !-- output rank

    !-- Pointers

    S % iCycle        =>  iCycle
    S % T             =>  T
    S % dT_Candidate  =>  dT_Candidate

    !-- Cleanup

    end associate !-- dT
    end associate !-- T_Max, etc.
    end associate !-- B

  end subroutine Initialize_B


  function TimerRecord ( S, Level ) result ( T )

    class ( Series_B_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimerRecord, &
               Name = trim ( S % Name ) // '_Rcd', &
               Level = Level )

  end function TimerRecord


  function TimerWrite ( S, Level ) result ( T )

    class ( Series_B_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimerWrite, &
               Name = trim ( S % Name ) // '_Wrt', &
               Level = Level )

  end function TimerWrite


  function TimerRead ( S, Level ) result ( T )

    class ( Series_B_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimerRead, &
               Name = trim ( S % Name ) // '_Rd', &
               Level = Level )

  end function TimerRead


  subroutine Record ( S )

    class ( Series_B_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iT  !-- iTimer

    associate &
      ( BV => S % Basic % Value, &
        T_Max   => S % TimerMax, &
        TV_Max  => S % TimerMax % Value, &
        TV_Min  => S % TimerMin % Value, &
        TV_Mean => S % TimerMean % Value, &
        dTV     => S % dT % Value, &
        iR      => S % iRecord )

    iR  =  iR + 1

    call Show ( 'Recording Series data', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )
    call Show ( iR, 'iRecord', S % IGNORABILITY )

    if ( iR  >  size ( BV, dim = 1 ) ) then
      call Show ( 'Too many time series entries', CONSOLE % ERROR )
      call Show ( 'Series_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Record', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    BV ( iR, S % TIME )   =  S % T
    BV ( iR, S % CYCLE )  =  S % iCycle

    associate ( MU  =>  PROGRAM_HEADER % MemoryUsage )
    BV ( iR, S % MEMORY_MAX_HWM )   =  MU % HighWaterMarkMax
    BV ( iR, S % MEMORY_MIN_HWM )   =  MU % HighWaterMarkMin
    BV ( iR, S % MEMORY_MEAN_HWM )  =  MU % HighWaterMarkMean
    BV ( iR, S % MEMORY_MAX_RSS )   =  MU % ResidentSetSizeMax
    BV ( iR, S % MEMORY_MIN_RSS )   =  MU % ResidentSetSizeMin
    BV ( iR, S % MEMORY_MEAN_RSS )  =  MU % ResidentSetSizeMean
    end associate !-- MU

    if ( S % iCycle  >  0 ) then

      associate &
        ( T_1D  =>  PROGRAM_HEADER % Timer_1D )
      associate &
        ( TimeMax   =>  T_1D % PreviousMax   +  T_1D % TimeMax, &
          TimeMin   =>  T_1D % PreviousMin   +  T_1D % TimeMin, &
          TimeMean  =>  T_1D % PreviousMean  +  T_1D % TimeMean )
      do iT  =  1, S % N_SERIES_TIMER
        TV_Max  ( iR, iT ) = TimeMax  ( iT ) / S % iCycle
        TV_Min  ( iR, iT ) = TimeMin  ( iT ) / S % iCycle
        TV_Mean ( iR, iT ) = TimeMean ( iT ) / S % iCycle
        call Show ( TV_Max ( iR, iT ), &
                    T_Max % Unit ( iT ), &
                    trim ( T_Max % Variable ( iT ) ) // ' per cycle', &
                    S % IGNORABILITY + 1 )
      end do !-- iT
      end associate !-- TimeMax, etc.
      end associate !-- T_1D

      dTV ( iR, : )  =  S % dT_Candidate

    end if

    end associate !-- BV, etc.

  end subroutine Record


  subroutine Write ( S )

    class ( Series_B_Form ), intent ( inout ) :: &
      S

    if ( .not. allocated ( S % GridImageStream ) ) &
      return

    associate &
      ( CI   =>  S % CurveImage, &
        GIS  =>  S % GridImageStream, &
        B    =>  S % Basic )
        
    !-- Assumes GridImageStream to which CI is attached is open

    call GIS % Open ( GIS % ACCESS_CREATE, SeriesOption = .false. )
    call Show ( 'Writing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )
    call CI % ClearGrid ( )
    call CI % SetGridWrite  &
           ( Directory = 'Series', &
             NodeCoordinate = B % Value ( 1 : S % iRecord, S % TIME ), &
             nProperCells = S % iRecord, &
             oValue = 0, &
             CoordinateUnitOption = B % Unit ( S % TIME ), &
             CoordinateLabelOption = 't' )
    call CI % Write ( )
    call GIS % Close ( ) 

    end associate !-- CI, etc.

  end subroutine Write


  subroutine Read ( S, nSeries )

    class ( Series_B_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      nSeries

    S % iRecord  =  nSeries

    if ( .not. allocated ( S % GridImageStream ) ) &
      return

    call Show ( 'Reading ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    associate &
      ( GIS => S % GridImageStream, &
        SB  => S % Basic, &
        CI => S % CurveImage )

    call GIS % Open ( GIS % ACCESS_READ, SeriesOption = .false. )
    call CI % ClearGrid ( )
    call CI % SetGridRead  &
           ( Directory = 'Series', &
             nProperCells = nSeries, &
             oValue = 0 )
    call CI % Read ( StorageOnlyOption = .true. )
    call GIS % Close ( ) 

    end associate !-- GIS, etc.

  end subroutine Read


  subroutine Restore ( S, iCycleRestart )

    class ( Series_B_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iCycleRestart

    integer ( KDI ) :: &
      iT  !-- iTimer

    associate &
      ( T_Mean   =>  S % TimerMean, &
        TV_Max   =>  S % TimerMax % Value, &
        TV_Min   =>  S % TimerMin % Value, &
        TV_Mean  =>  S % TimerMean % Value, &
        iR  =>  S % iRecord )

    call Show ( 'Restoring Series data', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )
    call Show ( iR, 'iRecord', S % IGNORABILITY )

    associate &
        ( TimeMax   =>  PROGRAM_HEADER % Timer_1D % PreviousMax, &
          TimeMin   =>  PROGRAM_HEADER % Timer_1D % PreviousMin, &
          TimeMean  =>  PROGRAM_HEADER % Timer_1D % PreviousMean )
    do iT  =  1,  S % N_SERIES_TIMER
      TimeMax  ( iT )  =  TV_Max  ( iR, iT )  *  iCycleRestart
      TimeMin  ( iT )  =  TV_Min  ( iR, iT )  *  iCycleRestart
      TimeMean ( iT )  =  TV_Mean ( iR, iT )  *  iCycleRestart
      call Show ( TimeMean ( iT ), &
                  T_Mean % Unit ( iT ), &
                  trim ( T_Mean % Variable ( iT ) ) // ' (TimeMean)', &
                  S % IGNORABILITY )
    end do !-- iT
    end associate !-- TimeMax, etc.

    end associate !-- T_Mean, etc.

  end subroutine Restore


  impure elemental subroutine Finalize ( S )

    type ( Series_B_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % CurveImage ) ) &
      deallocate ( S % CurveImage )
    if ( allocated ( S % GridImageStream ) ) &
      deallocate ( S % GridImageStream )
    if ( allocated ( S % dT ) ) &
      deallocate ( S % dT )
    if ( allocated ( S % TimerMean ) ) &
      deallocate ( S % TimerMean )
    if ( allocated ( S % TimerMin ) ) &
      deallocate ( S % TimerMin )
    if ( allocated ( S % TimerMax ) ) &
      deallocate ( S % TimerMax )
    if ( allocated ( S % Basic ) ) &
      deallocate ( S % Basic )

    if ( S % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

  end subroutine Finalize


end module Series_B__Form
