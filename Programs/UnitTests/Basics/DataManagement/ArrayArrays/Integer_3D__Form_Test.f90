program Integer_3D__Form_Test

  use Specifiers
  use Integer_3D__Form

  implicit none
  
  integer ( KDI ), dimension ( 3 ), parameter :: &
    SHAPE = [ 4, 3, 2 ]
  integer ( KDI ) :: &
    i
  type ( Integer_3D_Form ), dimension ( 3 ) :: &
    I_3D

  call I_3D ( 1 ) % Initialize ( SHAPE, ClearOption = .true. )
  print *
  print *, 'I_3D ( 1 ) % Value = ', I_3D ( 1 ) % Value 

  call I_3D ( 2 ) % Initialize ( SHAPE, iaLowerBoundOption = [ -1, -1, -1 ] )
  I_3D ( 2 ) % Value = reshape ( [ ( i, i = 1, product ( SHAPE ) ) ], SHAPE )
  print *
  print *, &
    'lbound ( I_3D ( 2 ) % Value ) = ', lbound ( I_3D ( 2 ) % Value )
  print *, 'I_3D ( 2 ) % Value = ', I_3D ( 2 ) % Value

  call I_3D ( 3 ) % Initialize ( I_3D ( 2 ) )
  print *
  print *, &
    'lbound ( I_3D ( 3 ) % Value ) = ', lbound ( I_3D ( 3 ) % Value )
  print *, 'I_3D ( 3 ) % Value = ', I_3D ( 3 ) % Value

end program Integer_3D__Form_Test
