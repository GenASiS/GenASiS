module Timer_1D__Form

  use Specifiers
  use Display
  use MessagePassing
  use Timer_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_TIMERS = 128

  type, public :: Timer_1D_Form
    integer ( KDI ) :: &
      MAX_TIMERS = MAX_TIMERS
    integer ( KDI ) :: &
      nTimers = 0, &
      LevelMin = 0, &
      LevelMax = 0
    real ( KDR ) :: &
      DisplayFraction
    real ( KDR ), dimension ( MAX_TIMERS ) :: &
      TimeThis = 0.0_KDR, &
      TimeMax  = 0.0_KDR, &
      TimeMin  = 0.0_KDR, &
      TimeMean = 0.0_KDR
    real ( KDR ), dimension ( MAX_TIMERS ) :: &
      PreviousMax  = 0.0_KDR, &  !-- From prior run segments before restart
      PreviousMin  = 0.0_KDR, &
      PreviousMean = 0.0_KDR
    type ( TimerForm ), dimension ( MAX_TIMERS ) :: &
      Element
    type ( TimerForm ), dimension ( MAX_TIMERS ) :: &
      TimerMax, &
      TimerMin, &
      TimerMean
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Get
    procedure, public, pass :: &
      Read
    final :: &
      Finalize
  end type Timer_1D_Form


contains


  subroutine Initialize ( T_1D, DisplayFraction, LevelMin )

    class ( Timer_1D_Form ), intent ( inout ) :: &
      T_1D
    real ( KDR ), intent ( in ) :: &
      DisplayFraction
    integer ( KDI ), intent ( in ) :: &
      LevelMin

    call Show ( 'Setting Timer parameters', CONSOLE % INFO_1 )

    T_1D % LevelMin  =  LevelMin
    call Show ( T_1D % LevelMin, 'LevelMin', CONSOLE % INFO_1 )

    T_1D % LevelMax  =  10  !-- Should match "Suffix" length in Timer_Form
    call Show ( T_1D % LevelMax, 'LevelMax', CONSOLE % INFO_1 )

    T_1D % DisplayFraction  =  DisplayFraction
    call Show ( T_1D % DisplayFraction, 'DisplayFraction', CONSOLE % INFO_1 )

  end subroutine Initialize


  subroutine Get ( T_1D, Handle, Name, Level, T )

    class ( Timer_1D_Form ), intent ( inout ), target :: &
      T_1D
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ) :: &
      Level    
    type ( TimerForm ), pointer, intent ( out ) :: &
      T

    T  =>  null ( )

    if ( Handle  <=  0 ) then

      associate ( nT  =>  T_1D % nTimers )
          nT  =  nT + 1
      Handle  =  nT
      end associate !-- nT

      associate ( T_New  =>  T_1D % Element ( Handle ) )
      call T_New % Initialize ( Name, Level, HandleOption = Handle )
      call Show ( 'Adding a Timer', CONSOLE % INFO_3 )
      call Show ( T_New % Name, 'Name', CONSOLE % INFO_3 )
      call Show ( T_New % Level, 'Level', CONSOLE % INFO_3 )
      call Show ( T_New % Handle, 'Handle', CONSOLE % INFO_3 )
      end associate !-- T_New

    end if

    T  =>  T_1D % Element ( Handle )

  end subroutine Get


  subroutine Read ( T_1D, Ignorability, CommunicatorOption )

    class ( Timer_1D_Form ), intent ( inout ) :: &
      T_1D
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
      
    integer ( KDI ) :: &
      iT
    real ( KDR ) :: &
      ExecutionTime, &
      MaxMinusMeanFraction
    type ( QuantityForm ) :: &
      MaxMinusMean
    logical ( KDL ) :: &
      Running
    type ( CollectiveOperation_R_Form ) :: &
      CO

    associate ( nT  =>  T_1D % nTimers )

    call Show ( 'This process timers', Ignorability + 1 )
    do iT  =  1,  nT
      associate &
        ( T         =>  T_1D % Element ( iT ), &
          TimeThis  =>  T_1D % TimeThis ( iT ) )
      Running  =  .false.
      if ( T % iStart  >  0 ) then
        Running  =  .true.
        call T % Stop ( )
      end if
      call T % ShowTotal ( Ignorability + 1 )
      TimeThis  =  T % TotalTime
      if ( Running ) &
        call T % Start ( )
      end associate !-- T
    end do !-- iT

    if ( present ( CommunicatorOption ) ) then

      call CO % Initialize &
             ( CommunicatorOption, nOutgoing = [ nT ], nIncoming = [ nT ] )

      CO % Outgoing % Value  =  T_1D % TimeThis ( 1 : nT )

      do iT  =  1,  nT
        associate &
          ( T       =>  T_1D % Element ( iT ), &
            T_Max   =>  T_1D % TimerMax ( iT ), &
            T_Min   =>  T_1D % TimerMin ( iT ), &
            T_Mean  =>  T_1D % TimerMean ( iT ) )
        call T_Max  % Initialize ( T )
        call T_Min  % Initialize ( T )
        call T_Mean % Initialize ( T )
        end associate !-- T, etc.
      end do !-- iT

      call Show ( 'Across processes max timers', Ignorability + 1 )
      call CO % Reduce ( REDUCTION % MAX )
      do iT  =  1,  nT
        associate &
          ( T_Max    =>  T_1D % TimerMax ( iT ), &
            TimeMax  =>  T_1D % TimeMax ( iT ) )
        call T_Max % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) )
        call T_Max % ShowTotal ( Ignorability + 1 )
        TimeMax  =  T_Max % TotalTime
        end associate !-- T_Max, etc.
      end do !-- iT
      
      call Show ( 'Across processes min timers', Ignorability + 1 )
      call CO % Reduce ( REDUCTION % MIN )
      do iT  =  1,  nT
        associate &
          ( T_Min    =>  T_1D % TimerMin ( iT ), &
            TimeMin  =>  T_1D % TimeMin ( iT ) )
        call T_Min % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) )
        call T_Min % ShowTotal ( Ignorability + 1 )
        TimeMin  =  T_Min % TotalTime
        end associate !-- T_Min, etc.
      end do !-- iT
      
      call Show ( 'Across processes mean timers', Ignorability + 1 )
      call CO % Reduce ( REDUCTION % SUM )
      do iT  =  1,  nT
        associate &
          ( T_Mean    =>  T_1D % TimerMean ( iT ), &
            TimeMean  =>  T_1D % TimeMean ( iT ) )
        call T_Mean % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT )  &
                      /  CommunicatorOption % Size )
        call T_Mean % ShowTotal ( Ignorability + 1 )
        TimeMean  =  T_Mean % TotalTime
        end associate !-- T_Mean, etc.
      end do !-- iT
      
      call Show ( 'Significant mean timers', Ignorability )
      do iT  =  1,  nT
        associate ( T_Mean  =>  T_1D % TimerMean ( iT ) )
        if ( iT  ==  1 ) &
          ExecutionTime  =  T_Mean % TotalTime
        if ( T_Mean % Level  <=  T_1D % LevelMin  &
             .or. ( T_Mean % TotalTime  /  ExecutionTime  &
                    >=  T_1D % DisplayFraction ) ) &
        then
          call T_Mean % ShowTotal ( Ignorability )
        end if
        end associate !-- T_Mean
      end do !-- iT

      call Show ( 'Significant MAX - MEAN', Ignorability )
      do iT  =  1, nT
        associate &
          ( T_Max   =>  T_1D % TimerMax ( iT ), &
            T_Mean  =>  T_1D % TimerMean ( iT ) )
        MaxMinusMean  =  T_Max % TotalTime  -  T_Mean % TotalTime
        if ( iT == 1 ) &
          ExecutionTime  = T_Mean % TotalTime
        if ( T_Mean % TotalTime  /  ExecutionTime  &
             >=  T_1D % DisplayFraction ) &
        then
          MaxMinusMeanFraction  &
            =  abs ( MaxMinusMean % Number ) &
                    /  max ( T_Mean % TotalTime % Number, &
                             sqrt ( tiny ( 0.0_KDR ) ) )
          if ( MaxMinusMeanFraction  >=  T_1D % DisplayFraction ) &
            call Show ( MaxMinusMean, &
                        'MAX - MEAN ' // T_Mean % Name, &
                        Ignorability )
        end if
        end associate !-- T_Max, etc.
      end do !-- iT

    end if !-- CommunicatorOption

    end associate !-- nT

  end subroutine Read


  impure elemental subroutine Finalize ( T_1D )

    type ( Timer_1D_Form ), intent ( inout ) :: &
      T_1D

  end subroutine Finalize

end module Timer_1D__Form
