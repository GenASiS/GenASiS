program Communicator_Form_Test

  use Display
  use Communicator_Form

  implicit none
  type ( CommunicatorForm ), allocatable :: &
    C, &
    SC

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  
  allocate ( SC )
  call SC % Initialize ( C, [ C % Rank ], NameOption = 'Subcommunicator' )
  
  deallocate ( SC )
  deallocate ( C )

end program Communicator_Form_Test
