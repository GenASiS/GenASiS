program Integer_2D__Form_Test

  use Specifiers
  use Integer_2D__Form

  implicit none
  
  integer ( KDI ), dimension ( 2 ), parameter :: &
    SHAPE = [ 4, 3 ]
  integer ( KDI ) :: &
    i
  type ( Integer_2D_Form ), dimension ( 3 ) :: &
    I_2D

  call I_2D ( 1 ) % Initialize ( SHAPE, ClearOption = .true. )
  print *
  print *, 'I_2D ( 1 ) % Value = ', I_2D ( 1 ) % Value 

  call I_2D ( 2 ) % Initialize ( SHAPE, iaLowerBoundOption = [ -1, -1 ] )
  I_2D ( 2 ) % Value = reshape ( [ ( i, i = 1, product ( SHAPE ) ) ], SHAPE )
  print *
  print *, &
    'lbound ( I_2D ( 2 ) % Value ) = ', lbound ( I_2D ( 2 ) % Value )
  print *, 'I_2D ( 2 ) % Value = ', I_2D ( 2 ) % Value

  call I_2D ( 3 ) % Initialize ( I_2D ( 2 ) )
  print *
  print *, &
    'lbound ( I_2D ( 3 ) % Value ) = ', lbound ( I_2D ( 3 ) % Value )
  print *, 'I_2D ( 3 ) % Value = ', I_2D ( 3 ) % Value

end program Integer_2D__Form_Test
