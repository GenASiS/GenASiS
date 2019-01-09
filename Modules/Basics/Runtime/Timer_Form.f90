module Timer_Form

  use Specifiers
  use Display
  use WallTime_Function

  implicit none
  private

  type, public :: TimerForm
    integer ( KDI ) :: &
      Level
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
    procedure, private, pass :: &
      InitializeNameLevel
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeNameLevel, InitializeClone
    procedure, public, pass :: &
      Start
    procedure, public, pass :: &
      Stop
    procedure, public, pass :: &
      ShowInterval
    procedure, public, pass :: &
      ShowTotal
  end type TimerForm

    character ( 8 ), private, parameter :: &
      Suffix = '::::::::'
    
contains


  subroutine InitializeNameLevel ( T, Name, Level )

    class ( TimerForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ) :: &
      Level

    T % Name  = Name
    T % Level = Level

    call T % StartTime % Initialize ( 's', 0.0_KDR )
    call T % StopTime % Initialize ( 's', 0.0_KDR )
    call T % TimeInterval % Initialize ( 's', 0.0_KDR )
    call T % TotalTime % Initialize ( 's', 0.0_KDR )

  end subroutine InitializeNameLevel


  subroutine InitializeClone ( T, T_Target )

    class ( TimerForm ), intent ( inout ) :: &
      T
    class ( TimerForm ), intent ( in ) :: &
      T_Target

    T % Name  = T_Target % Name
    T % Level = T_Target % Level

    call T % StartTime % Initialize ( 's', 0.0_KDR )
    call T % StopTime % Initialize ( 's', 0.0_KDR )
    call T % TimeInterval % Initialize ( 's', 0.0_KDR )
    call T % TotalTime % Initialize ( 's', 0.0_KDR )

  end subroutine InitializeClone


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

    call Show ( T % TimeInterval, &
                trim ( T % Name ) // ' Interval ' // Suffix ( 1 : T % Level ), &
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

    call Show ( T % TotalTime, &
                trim ( T % Name ) // ' ' // Suffix ( 1 : T % Level ), &
                Ignorability )

  end subroutine ShowTotal


end module Timer_Form
