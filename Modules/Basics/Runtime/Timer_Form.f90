module Timer_Form

  use VariableManagement
  use Display
  use WallTime_Function

  implicit none
  private

  type, public :: TimerForm
    type ( MeasuredValueForm ) :: &
      StartTime, &
      StopTime, &
      TimeInterval, &
      TotalTime
    logical ( KDL ) :: &
      Running = .false.
    character ( LDL ) :: &
      Name = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Start
    procedure, public, pass :: &
      Stop
    procedure, public, pass :: &
      ShowInterval
    procedure, public, pass :: &
      ShowTotal
  end type TimerForm

contains


  subroutine Initialize ( T, Name )

    class ( TimerForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in ) :: &
      Name

    T % Name = Name

    call T % StartTime % Initialize ( 's', 0.0_KDR )
    call T % StopTime % Initialize ( 's', 0.0_KDR )
    call T % TimeInterval % Initialize ( 's', 0.0_KDR )
    call T % TotalTime % Initialize ( 's', 0.0_KDR )

  end subroutine Initialize


  subroutine Start ( T )

    class ( TimerForm ), intent ( inout ) :: &
      T

    if ( T % Running ) return

    T % StartTime = WallTime ( )

    T % Running = .true.

  end subroutine Start


  subroutine Stop ( T )

    class ( TimerForm ), intent ( inout ) :: &
      T

    if ( .not. T % Running ) return

    T % StopTime = WallTime ( )

    T % TimeInterval  =  T % StopTime   -  T % StartTime
    T % TotalTime     =  T % TotalTime  +  T % TimeInterval

    T % Running = .false.

  end subroutine Stop


  impure elemental subroutine ShowInterval ( T, IgnorabilityOption )

    class ( TimerForm ), intent ( inout ) :: &
      T
    integer, intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      Ignorability

    if ( T % Name == '' ) return

    Ignorability = CONSOLE % INFO_2
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    call Show ( T % TimeInterval, trim ( T % Name ) // ' TimeInterval', &
                Ignorability )

  end subroutine ShowInterval


  impure elemental subroutine ShowTotal ( T, IgnorabilityOption )

    class ( TimerForm ), intent ( inout ) :: &
      T
    integer, intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      Ignorability

    if ( T % Name == '' ) return

    Ignorability = CONSOLE % INFO_2
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    call Show ( T % TotalTime, trim ( T % Name ) // ' TotalTime', &
                Ignorability )

  end subroutine ShowTotal


end module Timer_Form
