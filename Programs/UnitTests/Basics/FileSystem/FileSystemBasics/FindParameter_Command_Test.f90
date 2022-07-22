program FindParameter_Command_Test

  use ISO_FORTRAN_ENV
  use Specifiers
  use Display
  use MessagePassing
  use FindParameter_Command

  implicit none

  integer ( KDI ) :: &
    Integer_0D
  integer ( KDI ), dimension ( 10 ) :: &
    Integer_1D
  real ( KDR ) :: &
    Real_0D, &
    RealUnit_0D
  real ( KDR ), dimension ( 10 ) :: &
    Real_1D, &
    RealUnit_1D
  type ( QuantityForm ) :: &
    Quantity_0D
  type ( QuantityForm ), dimension ( 10 ) :: &
    Quantity_1D
  logical ( KDL ) :: &
    Logical_0D
  logical ( KDL ), dimension ( 10 ) :: &
    Logical_1D
  character ( 5 ) :: &
    Encoding
  character ( LDF ) :: &
    Character_0D
  character ( LDF ), dimension ( 10 ) :: &
    Character_1D
  character ( LDF ), dimension ( 13 ) :: &
    Buffer 
  type ( CommunicatorForm ), allocatable :: &
    C

!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_4' )

  call UNIT % Initialize ( )

  Buffer ( 1 ) = '!-- Commented string should not be read'
  Buffer ( 2 ) = 'ScalarInteger=10'
  Buffer ( 3 ) = 'ScalarReal=10.1'
  Buffer ( 4 ) = 'ScalarRealUnit=10.0~KILOMETER'
  Buffer ( 5 ) = 'ScalarQuantity=20.0~KILOMETER'
  Buffer ( 6 ) = 'ScalarLogical=True'
  Buffer ( 7 ) = 'ScalarString=Hello World'
  Buffer ( 8 ) = 'ArrayInteger=10,20,30'
  Buffer ( 9 ) = 'ArrayReal=10.5,20.5,30.5,40.5,50.5'
  Buffer ( 10 )  &
    = 'ArrayRealUnit=10.5~SECOND,20.5,30.5~MILLISECOND,40.5,50.5~MILLISECOND'
  Buffer ( 11 )  &
    = 'ArrayQuantity=10.5~SECOND,20.5,30.5~MILLISECOND,40.5,' &
      // '50.5~MILLISECOND'
  Buffer ( 12 ) = 'ArrayLogical=T,F,False,True,True,F'
  Buffer ( 13 ) = 'ArrayString=Lorem, ipsum, dolor, sit, amet'

  call FindParameter ( Integer_0D, Buffer, 'Test Buffer', 'ScalarInteger' )
  call FindParameter ( Real_0D, Buffer, 'TestBuffer', 'ScalarReal' )
  call FindParameter ( RealUnit_0D, Buffer, 'TestBuffer', 'ScalarRealUnit' )
  call FindParameter &
         ( Quantity_0D, Buffer, 'TestBuffer', 'ScalarQuantity' )
  call FindParameter ( Logical_0D, Buffer, 'Test Buffer', 'ScalarLogical' )
  call FindParameter ( Character_0D, Buffer, 'Test Buffer', 'ScalarString' )

  call FindParameter ( Integer_1D, Buffer, 'Test Buffer', 'ArrayInteger' )
  call FindParameter ( Real_1D, Buffer, 'TestBuffer', 'ArrayReal' )
  call FindParameter ( RealUnit_1D, Buffer, 'TestBuffer', 'ArrayRealUnit' )
  call FindParameter &
         ( Quantity_1D, Buffer, 'TestBuffer', 'ArrayQuantity' )
  call FindParameter ( Logical_1D, Buffer, 'Test Buffer', 'ArrayLogical' )
  call FindParameter ( Character_1D, Buffer, 'Test Buffer', 'ArrayString' )

  deallocate ( C )

end program FindParameter_Command_Test
