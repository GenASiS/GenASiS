program Clear_Command_Test

  use Specifiers
  use Clear_Command

  implicit none

  integer ( KDI ), dimension ( 2 ) :: &
    Integer_1D
  integer ( KDI ), dimension ( 2, 2 ) :: &
    Integer_2D
  integer ( KDI ), dimension ( 2, 2, 2 ) :: &
    Integer_3D
  integer ( KDI ), dimension ( 2, 2, 2, 2 ) :: &
    Integer_4D
  integer ( KBI ), dimension ( 2 ) :: &
    BigInteger_1D
  integer ( KBI ), dimension ( 2, 2 ) :: &
    BigInteger_2D
  integer ( KBI ), dimension ( 2, 2, 2 ) :: &
    BigInteger_3D
  integer ( KBI ), dimension ( 2, 2, 2, 2 ) :: &
    BigInteger_4D
  real ( KDR ), dimension ( 2 ) :: &
    Real_1D
  real ( KDR ), dimension ( 2, 2 ) :: &
    Real_2D
  real ( KDR ), dimension ( 2, 2, 2 ) :: &
    Real_3D
  real ( KDR ), dimension ( 2, 2, 2, 2 ) :: &
    Real_4D
  complex ( KDC ), dimension ( 2 ) :: &
    Complex_1D
  complex ( KDC ), dimension ( 2, 2 ) :: &
    Complex_2D
  complex ( KDC ), dimension ( 2, 2, 2 ) :: &
    Complex_3D
  complex ( KDC ), dimension ( 2, 2, 2, 2 ) :: &
    Complex_4D
  logical ( KDL ), dimension ( 2 ) :: &
    Logical_1D
  logical ( KDL ), dimension ( 2, 2 ) :: &
    Logical_2D
  logical ( KDL ), dimension ( 2, 2, 2 ) :: &
    Logical_3D
  logical ( KDL ), dimension ( 2, 2, 2, 2 ) :: &
    Logical_4D

  !-- Integer

  Integer_1D = 1_KDI
  Integer_2D = 1_KDI
  Integer_3D = 1_KDI
  Integer_4D = 1_KDI

  print *, 'Integer_1D = ', Integer_1D
  print *, 'Integer_2D = ', Integer_2D
  print *, 'Integer_3D = ', Integer_3D
  print *, 'Integer_4D = ', Integer_4D
  print *

  call Clear ( Integer_1D )
  call Clear ( Integer_2D )
  call Clear ( Integer_3D )
!  call Clear ( Integer_4D )

  print *, 'Integer_1D = ', Integer_1D
  print *, 'Integer_2D = ', Integer_2D
  print *, 'Integer_3D = ', Integer_3D
  print *, 'Integer_4D = ', Integer_4D
  print *

  !-- BigInteger

  BigInteger_1D = 1_KBI
  BigInteger_2D = 1_KBI
  BigInteger_3D = 1_KBI
  BigInteger_4D = 1_KBI

  print *, 'BigInteger_1D = ', BigInteger_1D
  print *, 'BigInteger_2D = ', BigInteger_2D
  print *, 'BigInteger_3D = ', BigInteger_3D
  print *, 'BigInteger_4D = ', BigInteger_4D
  print *

!  call Clear ( BigInteger_1D )
!  call Clear ( BigInteger_2D )
!  call Clear ( BigInteger_3D )
!  call Clear ( BigInteger_4D )

  print *, 'BigInteger_1D = ', BigInteger_1D
  print *, 'BigInteger_2D = ', BigInteger_2D
  print *, 'BigInteger_3D = ', BigInteger_3D
  print *, 'BigInteger_4D = ', BigInteger_4D
  print *

  !-- Real

  Real_1D = 1.0_KDR
  Real_2D = 1.0_KDR
  Real_3D = 1.0_KDR
  Real_4D = 1.0_KDR

  print *, 'Real_1D = ', Real_1D
  print *, 'Real_2D = ', Real_2D
  print *, 'Real_3D = ', Real_3D
  print *, 'Real_4D = ', Real_4D
  print *

  call Clear ( Real_1D )
  call Clear ( Real_2D )
  call Clear ( Real_3D )
!  call Clear ( Real_4D )

  print *, 'Real_1D = ', Real_1D
  print *, 'Real_2D = ', Real_2D
  print *, 'Real_3D = ', Real_3D
  print *, 'Real_4D = ', Real_4D
  print *

  !-- Complex

  Complex_1D = ( 1.0_KDR, 1.0_KDR )
  Complex_2D = ( 1.0_KDR, 1.0_KDR )
  Complex_3D = ( 1.0_KDR, 1.0_KDR )
  Complex_4D = ( 1.0_KDR, 1.0_KDR )

  print *, 'Complex_1D = ', Complex_1D
  print *, 'Complex_2D = ', Complex_2D
  print *, 'Complex_3D = ', Complex_3D
  print *, 'Complex_4D = ', Complex_4D
  print *

!  call Clear ( Complex_1D )
!  call Clear ( Complex_2D )
  call Clear ( Complex_3D )
!  call Clear ( Complex_4D )

  print *, 'Complex_1D = ', Complex_1D
  print *, 'Complex_2D = ', Complex_2D
  print *, 'Complex_3D = ', Complex_3D
  print *, 'Complex_4D = ', Complex_4D
  print *

  !-- Logical

  Logical_1D = .true.
  Logical_2D = .true.
  Logical_3D = .true.
  Logical_4D = .true.

  print *, 'Logical_1D = ', Logical_1D
  print *, 'Logical_2D = ', Logical_2D
  print *, 'Logical_3D = ', Logical_3D
  print *, 'Logical_4D = ', Logical_4D
  print *

  call Clear ( Logical_1D )
!  call Clear ( Logical_2D )
!  call Clear ( Logical_3D )
!  call Clear ( Logical_4D )

  print *, 'Logical_1D = ', Logical_1D
  print *, 'Logical_2D = ', Logical_2D
  print *, 'Logical_3D = ', Logical_3D
  print *, 'Logical_4D = ', Logical_4D
  print *

end program Clear_Command_Test
