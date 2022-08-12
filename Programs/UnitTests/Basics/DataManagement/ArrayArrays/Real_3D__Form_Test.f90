program Real_3D__Form_Test

  use Specifiers
  use Real_3D__Form

  implicit none
  
  integer ( KDI ), dimension ( 3 ), parameter :: &
    SHAPE = [ 4, 3, 2 ]
  integer ( KDI ) :: &
    i
  type ( Real_3D_Form ), dimension ( 5 ) :: &
    R_3D

  call R_3D ( 1 ) % Initialize ( SHAPE, ClearOption = .true. )
  print *
  print *, 'R_3D ( 1 ) % Value = ', R_3D ( 1 ) % Value 

  call R_3D ( 2 ) % Initialize ( SHAPE, iaLowerBoundOption = [ -1, -1, -1 ] )
  associate ( PI => CONSTANT % PI )
  R_3D ( 2 ) % Value &
    = reshape ( [ ( i * PI, i = 1, product ( SHAPE ) ) ], SHAPE )
  end associate
  print *
  print *, &
    'lbound ( R_3D ( 2 ) % Value ) = ', lbound ( R_3D ( 2 ) % Value )
  print *, 'R_3D ( 2 ) % Value = ', R_3D ( 2 ) % Value

  call R_3D ( 3 ) % Initialize ( R_3D ( 2 ) )
  print *
  print *, &
    'lbound ( R_3D ( 3 ) % Value ) = ', lbound ( R_3D ( 3 ) % Value )
  print *, 'R_3D ( 3 ) % Value = ', R_3D ( 3 ) % Value
  
  print*, 'R_3D ( 4 )'
  call R_3D ( 4 ) % Initialize ( [ 32, 64, 128 ], ClearOption = .true. )
  call random_number ( R_3D ( 4 ) % Value )
  call R_3D ( 4 ) % AllocateDevice ( )
  call R_3D ( 4 ) % UpdateDevice ( )
  print*, 'End R_3D ( 4 )'
  
  print*, 'R_3D ( 5 )'
  call R_3D ( 5 ) % Initialize ( [ 32, 64, 128 ], ClearOption = .true. )
  call random_number ( R_3D ( 5 ) % Value )
  call R_3D ( 5 ) % AllocateDevice ( )
  call R_3D ( 5 ) % UpdateDevice ( )
  print*, 'End R_3D ( 5 )'
  
  call R_3D ( 5 ) % MultiplyAdd ( R_3D ( 4 ), 13.0_KDR, &
                                  UseDeviceOption = .true. )

end program Real_3D__Form_Test
