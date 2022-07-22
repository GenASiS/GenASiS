program ParameterStream_Form_Test

  use Specifiers
  use Display
  use MessagePassing
  use ParameterStream_Form

  implicit none

  integer ( KDI ) :: &
    Integer
  integer ( KDI ), dimension ( 10 ) :: &
    IntegerArray_1, &
    IntegerArray_2, &
    IntegerArray_3
  real ( KDR ) :: &
    Real, &
    RealUnit
  real ( KDR ), dimension ( 10 ) :: &
    RealArray_1, &
    RealArray_2, &
    RealArray_3, &
    RealUnitArray_1, &
    RealUnitArray_2, &
    RealUnitArray_3
  type ( QuantityForm ) :: &
    Quantity
  type ( QuantityForm ), dimension ( 10 ) :: &
    QuantityArray_1, &
    QuantityArray_2, &
    QuantityArray_3
  logical ( KDL ) :: &
    Logical
  logical ( KDL ), dimension ( 10 ) :: & 
    LogicalArray_1, &
    LogicalArray_2, &
    LogicalArray_3
  character ( LDL ) :: & 
    String
  character ( LDL ), dimension ( 10 ) :: & 
    StringArray_1 = '', &
    StringArray_2 = '', &
    StringArray_3 = ''
  character ( LDF ) :: &
   Filename = 'ParameterStream_Form_Test_Parameters'
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( ParameterStreamForm ) :: &
    PS

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )

  call UNIT % Initialize ( )

  call PS % Initialize ( Filename, C % Rank )

  call PS % Read ( Integer, 'Integer' )
  call PS % Read ( IntegerArray_1, 'IntegerArray_1' )
  call PS % Read ( IntegerArray_2, 'IntegerArray_2' )
  call PS % Read ( IntegerArray_3, 'IntegerArray_3' )

  call PS % Read ( Real, 'Real' )
  call PS % Read ( RealArray_1, 'RealArray_1' )
  call PS % Read ( RealArray_2, 'RealArray_2' )
  call PS % Read ( RealArray_3, 'RealArray_3' )
  
  call PS % Read ( RealUnit, 'RealUnit' )
  call PS % Read ( RealUnitArray_1, 'RealUnitArray_1' )
  call PS % Read ( RealUnitArray_2, 'RealUnitArray_2' )
  call PS % Read ( RealUnitArray_3, 'RealUnitArray_3' )
  
  call PS % Read ( Quantity, 'Quantity' )
  call PS % Read ( QuantityArray_1, 'QuantityArray_1' )
  call PS % Read ( QuantityArray_2, 'QuantityArray_2' )
  call PS % Read ( QuantityArray_3, 'QuantityArray_3' )
  
  call PS % Read ( Logical, 'Logical' )
  call PS % Read ( LogicalArray_1, 'LogicalArray_1' )
  call PS % Read ( LogicalArray_2, 'LogicalArray_2' )
  call PS % Read ( LogicalArray_3, 'LogicalArray_3' )
  
  call PS % Read ( String, 'String' )
  call PS % Read ( StringArray_1, 'StringArray_1' )
  call PS % Read ( StringArray_2, 'StringArray_2' )
  call PS % Read ( StringArray_3, 'StringArray_3' )
  
  !--  Test for retaining default value
  call PS % Read ( Real, 'BogusReal' )
  call Show ( Real, 'Previously read real' ) 
  
  call PS % Read ( StringArray_3, 'MyRealArray' )
  call Show ( StringArray_3, 'MyRealArray default' )  

  deallocate ( C )
  
end program ParameterStream_Form_Test
