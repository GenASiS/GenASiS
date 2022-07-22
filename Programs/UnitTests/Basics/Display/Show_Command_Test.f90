program Show_Command_Test

  use iso_fortran_env
  use iso_c_binding
  use Specifiers
  use CONSOLE_Singleton
  use Show_Command

  implicit none

  include 'mpif.h'

  integer ( KDI ) :: &
    iV, &
    Rank, &
    Error
  character ( 5 ) :: &
    Encoding
  type ( QuantityForm ) :: &
    Length
  real ( KDR ), dimension ( 5 ), target :: &
    A
  type ( c_ptr ), dimension ( 5 ) :: &
    Addresses

!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  call MPI_INIT ( Error )
  call MPI_COMM_RANK ( MPI_COMM_WORLD, Rank, Error )
  
  call UNIT % Initialize ( )

  call CONSOLE % Initialize ( ProcessRank = Rank )
  call CONSOLE % SetVerbosity ( 'INFO_5' )
  call CONSOLE % SetDisplayRank ( 0 )

  call MPI_BARRIER ( MPI_COMM_WORLD, Error )

  call Show ( 'Tests of Show_Command' )
  call Show ( Rank, 'Rank' )
  
  call Show &
         ( 1_KDI, 'Test integer' )
  call Show &
         ( [ 2_KDI, 3_KDI, 4_KDI ], 'Test 1D integer array' )
  call Show &
         ( reshape ( [ 5_KDI, 6_KDI, 7_KDI, 8_KDI ], [ 2, 2 ] ), &
           'Test 2D integer array' )

  call Show &
         ( 10.0_KDR, 'Test real' )
  call Show &
         ( [ 20.0_KDR, 30.0_KDR, 40.0_KDR ], 'Test 1D real array' )
  call Show &
         ( reshape ( [ 50.0_KDR, 60.0_KDR, 70.0_KDR, 80.0_KDR ], [ 2, 2 ] ), &
           'Test 2D real array' )

  call Show ( .true., 'Test logical' )
  call Show ( [ .true., .true., .false. ], 'Test logical array' ) 

  call Show ( 'Hello world', 'Test character' )
  call Show &
         ( [ 'Hello world 1', 'Hello world 2', 'Hello world 3' ], &
           'Test character array' )

  Length = 10.0_KDR  *  UNIT % METER
  call Show ( Length, 'Length in program units' )
  call Show ( Length, UNIT % METER, 'Length in m' )  
  call Show ( Length, UNIT % KILOMETER, 'Length in km' )
  call Show ( Length, UNIT % SECOND,  'Length in second' )
  call Show ( spread ( Length % Number, 1, 10 ), UNIT % PARSEC, &
              'Length in pc' )
  
  Addresses = [ ( c_loc ( A ( iV ) ), iV = 1, size ( A ) ) ]
  call SHow ( Addresses, 'Elements of A' )
  !-- FIXME: This breaks on GCC 11.3.0. For some reason it does not resolve
  !          to Show_C_Pointer, but to Show_C_Pointer_1D, where it causes
  !          a seg fault.
  call Show ( c_loc ( A ), 'Address of A' )
  
  call MPI_FINALIZE ( Error )

end program Show_Command_Test
