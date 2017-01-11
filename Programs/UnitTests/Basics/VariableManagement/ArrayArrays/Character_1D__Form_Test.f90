program Character_1D__Form_Test

  use Specifiers
  use Character_1D__Form

  implicit none

  integer ( KDI ), parameter :: &
    SIZE = 4
  integer ( KDI ) :: &
    i
  type ( Character_1D_Form ) :: &
    C

  call C % Initialize ( SIZE, iLowerBoundOption = -1 )
  C % Value ( -1 ) = "These are the times"
  C % Value ( 0 ) = "Bad to the bone"
  C % Value ( 1 ) = "Mutha f***a!"
  C % Value ( 2 ) = "Are you kidding me?"
  print *
  print *, &
    'lbound ( C % Value ) = ', lbound ( C % Value )
  print *, 'C % Value = ', C % Value

end program Character_1D__Form_Test
