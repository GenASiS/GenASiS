program DelayFileAccess_Command_Test

  use Display
  use MessagePassing
  use DelayFileAccess_Command

  implicit none
  type ( CommunicatorForm ), allocatable :: &
    C

  allocate ( C )

  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )

  call DelayFileAccess ( C % Rank )

  deallocate ( C )

end program DelayFileAccess_Command_Test
