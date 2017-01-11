program CommandLineOptions_Form_Test

  use VariableManagement
  use Display
  use MessagePassing
  use CommandLineOptions_Form

  implicit none

  
  integer ( KDI ) :: &
    TestInteger
  integer ( KDI ), dimension ( 3 ) :: &
    TestIntegerArray
  real ( KDR ) :: &
    TestReal
  logical ( KDL ) :: &
    TestLogical
  character ( LDL ) :: &
    TestString, &
    OptionStringLabel  !-- FIXME: Needed for Intel 12.1.3
  character ( LDL ), dimension ( 3 ) :: &
    TestStringArray
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( CommandLineOptionsForm ) :: &
    CLO

  !-- Run this program using the following arguments:
  ! ./CommandLineOptions_Form_Test_Jaguar_Cray OptionInteger=2 \
  !     OptionReal=3.14 OptionIntegerArray=64,64,64 OptionLogical=true \
  !     OptionString=helloworld

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
    
  call UNIT % Initialize ( )

  call CLO % Initialize ( ) 

  TestInteger = 2
  call CLO % Read ( TestInteger, 'OptionInteger' )
  
  TestIntegerArray = 1
  call CLO % Read ( TestIntegerArray, 'OptionIntegerArray' )
  
  TestReal = 1.5
  call CLO % Read ( TestReal, 'OptionReal' )
  
  TestLogical = .false.
  call CLO % Read ( TestLogical, 'OptionLogical' )
  
  TestString = 'Undefined'
!-- FIXME: Intel 12.1.3 garbled this explicit argument
!  call CLO % Read ( TestString, 'OptionString' )
  OptionStringLabel = 'OptionString'
  call CLO % Read ( TestString, OptionStringLabel )
  
  TestStringArray ( 1 ) = 'One'
  TestStringArray ( 2 ) = 'Two'
  TestStringArray ( 3 ) = 'Three'
  call CLO % Read ( TestStringArray, 'OptionStringArray' )
   
  call Show ( TestInteger, 'OptionInteger' )
  call Show ( TestIntegerArray,  'OptionIntegerArray'    )
  call Show ( TestReal,    'OptionReal'    )
  call Show ( TestLogical, 'OptionLogical' )
  call Show ( TestString,  'OptionString'  )
  call Show ( TestStringArray,  'OptionStringArray'  )

  deallocate ( C )

end program CommandLineOptions_Form_Test
