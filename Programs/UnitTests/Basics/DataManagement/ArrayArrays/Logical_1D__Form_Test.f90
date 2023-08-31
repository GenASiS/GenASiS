program Logical_1D__Form_Test

  use Specifiers
  use ArrayOperations
  use Logical_1D__Form

  implicit none

  integer ( KDI ), parameter :: &
    SIZE = 4
  integer ( KDL ) :: &
    i
  type ( Logical_1D_Form ), dimension ( : ), allocatable :: &
    L_1D

  allocate ( L_1D ( 3 ) )
  
  call L_1D ( 1 ) % Initialize ( SIZE, ClearOption = .true. )
  print *
  print *, 'L_1D ( 1 ) % Value = ', L_1D ( 1 ) % Value 

  call L_1D ( 2 ) % Initialize ( SIZE, iLowerBoundOption = -1 )
  L_1D ( 2 ) % Value = [ ( mod ( i, 2 ) == 0, i = 1, SIZE ) ]
  print *
  print *, &
    'lbound ( L_1D ( 2 ) % Value ) = ', lbound ( L_1D ( 2 ) % Value )
  print *, 'L_1D ( 2 ) % Value = ', L_1D ( 2 ) % Value

  call L_1D ( 3 ) % Initialize ( L_1D ( 2 ) )
  print *
  print *, &
    'lbound ( L_1D ( 3 ) % Value ) = ', lbound ( L_1D ( 3 ) % Value )
  print *, 'L_1D ( 3 ) % Value = ', L_1D ( 3 ) % Value
  
  call L_1D ( 2 ) % AllocateDevice ( )
  call L_1D ( 2 ) % UpdateDevice ( )
  
  call Clear ( L_1D ( 2 ) % Value )
  
  print *, &
    'lbound ( L_1D ( 2 ) % Value ) = ', lbound ( L_1D ( 2 ) % Value )
  print *, 'Cleared L_1D ( 2 ) % Value = ', L_1D ( 2 ) % Value
  
  call L_1D ( 2 ) % UpdateHost ( )
  print *, &
    'lbound ( L_1D ( 2 ) % Value ) = ', lbound ( L_1D ( 2 ) % Value )
  print *, 'Original L_1D ( 2 ) % Value = ', L_1D ( 2 ) % Value

  deallocate ( L_1D )

end program Logical_1D__Form_Test
