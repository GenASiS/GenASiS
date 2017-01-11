program WallTime_Function_Test

  use VariableManagement
  use Display
  use MessagePassing
  use WallTime_Function

  implicit none

  integer ( KDI ) :: &
    i
  type ( MeasuredValueForm ) :: &
    WallTimeStart, &
    WallTimeFinish, &
    WallTimeFinish2
  type ( CommunicatorForm ), allocatable :: &
    C

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )

  WallTimeStart = WallTime ( )

  do i = 1, 10000000
    !-- do nothing
  end do

  WallTimeFinish = WallTime ( )

  call Show ( WallTimeStart, 'WallTimeStart' )
  call Show ( WallTimeFinish, 'WallTimeFinish' )
  call Show ( WallTimeFinish - WallTimeStart, 'Elapsed time' )

  deallocate ( C )

end program WallTime_Function_Test
