program Real_2D__Form_Test

  use Specifiers
  use Real_2D__Form

  implicit none
  
  integer ( KDI ), dimension ( 2 ), parameter :: &
    SHAPE = [ 4, 3 ]
  integer ( KDI ) :: &
    i
  type ( Real_2D_Form ), dimension ( 3 ) :: &
    R_2D

  call R_2D ( 1 ) % Initialize ( SHAPE, ClearOption = .true. )
  print *
  print *, 'R_2D ( 1 ) % Value = ', R_2D ( 1 ) % Value 

  call R_2D ( 2 ) % Initialize ( SHAPE, iaLowerBoundOption = [ -1, -1 ] )
  associate ( PI => CONSTANT % PI )
  R_2D ( 2 ) % Value &
    = reshape ( [ ( i * PI, i = 1, product ( SHAPE ) ) ], SHAPE )
  end associate
  print *
  print *, &
    'lbound ( R_2D ( 2 ) % Value ) = ', lbound ( R_2D ( 2 ) % Value )
  print *, 'R_2D ( 2 ) % Value = ', R_2D ( 2 ) % Value

  call R_2D ( 3 ) % Initialize ( R_2D ( 2 ) )
  print *
  print *, &
    'lbound ( R_2D ( 3 ) % Value ) = ', lbound ( R_2D ( 3 ) % Value )
  print *, 'R_2D ( 3 ) % Value = ', R_2D ( 3 ) % Value

end program Real_2D__Form_Test
