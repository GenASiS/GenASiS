program Real_1D__Form_Test

  use Specifiers
  use Real_1D__Form

  implicit none

  integer ( KDI ), parameter :: &
    SIZE = 4
  integer ( KDI ) :: &
    i
  type ( Real_1D_Form ), dimension ( 3 ) :: &
    R_1D

  call R_1D ( 1 ) % Initialize ( SIZE, ClearOption = .true. )
  print *
  print *, 'R_1D ( 1 ) % Value = ', R_1D ( 1 ) % Value 

  call R_1D ( 2 ) % Initialize ( SIZE, iLowerBoundOption = -1 )
  associate ( PI => CONSTANT % PI )
  R_1D ( 2 ) % Value = [ ( i * PI, i = 1, SIZE ) ]
  end associate
  print *
  print *, &
    'lbound ( R_1D ( 2 ) % Value ) = ', lbound ( R_1D ( 2 ) % Value )
  print *, 'R_1D ( 2 ) % Value = ', R_1D ( 2 ) % Value

  call R_1D ( 3 ) % Initialize ( R_1D ( 2 ) )
  print *
  print *, &
    'lbound ( R_1D ( 3 ) % Value ) = ', lbound ( R_1D ( 3 ) % Value )
  print *, 'R_1D ( 3 ) % Value = ', R_1D ( 3 ) % Value

end program Real_1D__Form_Test
