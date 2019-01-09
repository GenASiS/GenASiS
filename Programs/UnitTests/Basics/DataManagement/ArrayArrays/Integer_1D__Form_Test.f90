program Integer_1D__Form_Test

  use Specifiers
  use Integer_1D__Form

  implicit none

  integer ( KDI ), parameter :: &
    SIZE = 4
  integer ( KDI ) :: &
    i
  type ( Integer_1D_Form ), dimension ( 3 ) :: &
    I_1D

  call I_1D ( 1 ) % Initialize ( SIZE, ClearOption = .true. )
  print *
  print *, 'I_1D ( 1 ) % Value = ', I_1D ( 1 ) % Value 

  call I_1D ( 2 ) % Initialize ( SIZE, iLowerBoundOption = -1 )
  I_1D ( 2 ) % Value = [ ( i, i = 1, SIZE ) ]
  print *
  print *, &
    'lbound ( I_1D ( 2 ) % Value ) = ', lbound ( I_1D ( 2 ) % Value )
  print *, 'I_1D ( 2 ) % Value = ', I_1D ( 2 ) % Value

  call I_1D ( 3 ) % Initialize ( I_1D ( 2 ) )
  print *
  print *, &
    'lbound ( I_1D ( 3 ) % Value ) = ', lbound ( I_1D ( 3 ) % Value )
  print *, 'I_1D ( 3 ) % Value = ', I_1D ( 3 ) % Value

end program Integer_1D__Form_Test
