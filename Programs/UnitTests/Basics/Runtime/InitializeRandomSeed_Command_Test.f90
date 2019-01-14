program InitializeRandomSeed_Command_Test

  use Specifiers
  use Display
  use MessagePassing
  use InitializeRandomSeed_Command
  
  implicit none
  
  type ( CommunicatorForm ), allocatable :: &
    C

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_4' )
  
  !-- dummy allocation to eat some memory
  
  call InitializeRandomSeed ( C )
  
  deallocate ( C )

end program InitializeRandomSeed_Command_Test
