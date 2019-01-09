program TableStream_Form_Test

  use Specifiers
  use Display
  use MessagePassing
  use TableStream_Form

  implicit none
  
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Value
  character ( LDF ) :: &
    Path     = '../Parameters/', &
    Filename = 'TableStream_Form_Test_Table'
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( TableStreamForm ) :: &
    TS

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_7' )

  call TS % Initialize ( Filename, C % Rank, PathOption = Path )
  
  call TS % Read ( Value )

  call Show ( Value, 'Table Values' )

  deallocate ( C )

end program TableStream_Form_Test
