program Timer_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use Timer_Form

  implicit none

  integer ( KDI ) :: &
    i
  type ( TimerForm ), dimension ( 5 ) :: &
    Timer
  type ( CommunicatorForm ), allocatable :: &
    C

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )

  call Timer ( 1 ) % Initialize ( 'TestTimer 1', Level = 1 )
  call Timer ( 2 ) % Initialize ( 'TestTimer 2', Level = 1 )

  !-- Interval A1; Timer ( 1 )

  call Show ( 'Interval A1', CONSOLE % INFO_3 )
  call Timer ( 1 ) % Start ( )

  do i = 1, 10000000
    !-- do nothing
  end do

  call Timer ( 1 ) % Stop ( )
  call Timer ( 1 ) % ShowInterval ( CONSOLE % INFO_3 )
  call Timer ( 1 ) % ShowTotal ( CONSOLE % INFO_3 )

  !-- Interval A2; Timer ( 1 )

  call Show ( 'Interval A2', CONSOLE % INFO_3 )
  call Timer ( 1 ) % Start ( )

  do i = 1, 20000000
    !-- do nothing
  end do

  call Timer ( 1 ) % Stop ( )
  call Timer ( 1 ) % ShowInterval ( CONSOLE % INFO_3 )
  call Timer ( 1 ) % ShowTotal ( CONSOLE % INFO_3 )

  !-- Interval B; Timer ( 2 )

  call Show ( 'Interval B', CONSOLE % INFO_3 )
  call Timer ( 2 ) % Start ( )

  do i = 1, 30000000
    !-- do nothing
  end do

  call Timer ( 2 ) % Stop ( )
  call Timer ( 2 ) % ShowInterval ( CONSOLE % INFO_3 )
  call Timer ( 2 ) % ShowTotal ( CONSOLE % INFO_3 )

  !-- Show Totals

  call Show ( 'Totals for all Timers' )
  call Timer % ShowTotal ( )

  deallocate ( C )

end program Timer_Form_Test
