program Copy_Command_Test

  use iso_c_binding
  use Specifiers
  use Devices
  use Clear_Command
  use Copy_Command

  implicit none

  integer ( KDI ), dimension ( 3 ) :: &
    Integer_1D, &
    Integer_1D_Copy
  integer ( KDI ), dimension ( 3, 3 ) :: &
    Integer_2D, &
    Integer_2D_Copy
  integer ( KDI ), dimension ( 3, 3, 3 ) :: &
    Integer_3D, &
    Integer_3D_Copy
  integer ( KBI ), dimension ( 3 ) :: &
    BigInteger_1D, &
    BigInteger_1D_Copy
  integer ( KBI ), dimension ( 3, 3 ) :: &
    BigInteger_2D, &
    BigInteger_2D_Copy
  integer ( KBI ), dimension ( 3, 3, 3 ) :: &
    BigInteger_3D, &
    BigInteger_3D_Copy
  real ( KDR ), dimension ( 3 ) :: &
    Real_1D, &
    Real_1D_Copy
  real ( KDR ), dimension ( 3, 3 ) :: &
    Real_2D, &
    Real_2D_Copy
  real ( KDR ), dimension ( 3, 3, 3 ) :: &
    Real_3D, &
    Real_3D_Copy
  complex ( KDC ), dimension ( 3 ) :: &
    Complex_1D, &
    Complex_1D_Copy
  complex ( KDC ), dimension ( 3, 3 ) :: &
    Complex_2D, &
    Complex_2D_Copy
  complex ( KDC ), dimension ( 3, 3, 3 ) :: &
    Complex_3D, &
    Complex_3D_Copy
  character ( 3 ) :: &
    Character_3
  character ( c_char ), dimension ( 4 ) :: &
    Character_1D
  type ( c_ptr ) :: &
    D_V, &
    D_V_Copy

  !-- Integer

  Integer_1D = 1_KDI
  call Clear ( Integer_1D_Copy )
  call Copy ( Integer_1D, Integer_1D_Copy )
  if ( all ( Integer_1D == Integer_1D_Copy ) ) then
    print*, 'Copy Integer_1D: Success' 
  else 
    print*, 'Copy Integer_1D: Fail' 
    print *, 'Integer_1D = ', Integer_1D
    print *, 'Integer_1D_Copy = ', Integer_1D_Copy
  end if
  print *
  print *

!  Integer_2D = 1_KDI
!  call Clear ( Integer_2D_Copy )
!  call Copy ( Integer_2D, Integer_2D_Copy )
!  print *, 'Integer_2D = ', Integer_2D
!  print *, 'Integer_2D_Copy = ', Integer_2D_Copy
!  print *
!  print *

!  Integer_3D = 1_KDI
!  call Clear ( Integer_3D_Copy )
!  print *, 'Integer_3D = ', Integer_3D
!  print *, 'Integer_3D_Copy = ', Integer_3D_Copy
!  print *

!  call Copy ( Integer_3D, Integer_3D_Copy )
!  print *, 'Integer_3D = ', Integer_3D
!  print *, 'Integer_3D_Copy = ', Integer_3D_Copy
!  print *
!  print *

!  !-- Integer section
!
!  Integer_1D = 1_KDI
!  call Clear ( Integer_1D_Copy )
!  print *, 'Integer_1D = ', Integer_1D
!  print *, 'Integer_1D_Copy = ', Integer_1D_Copy
!  print *
!
!  call Copy ( Integer_1D, 1, 1, 2, Integer_1D_Copy )
!  print *, 'Integer_1D = ', Integer_1D
!  print *, 'Integer_1D_Copy (Section ) = ', Integer_1D_Copy
!  print *
!  print *
!
!  Integer_2D = 1_KDI
!  call Clear ( Integer_2D_Copy )
!  print *, 'Integer_2D = ', Integer_2D
!  print *, 'Integer_2D_Copy = ', Integer_2D_Copy
!  print *
!
!  call Copy ( Integer_2D, [ 1, 1 ], [ 1, 1 ], [ 2, 2 ], Integer_2D_Copy )
!  print *, 'Integer_2D = ', Integer_2D
!  print *, 'Integer_2D_Copy (Section) = ', Integer_2D_Copy
!  print *
!  print *
!
!  Integer_3D = 1_KDI
!  call Clear ( Integer_3D_Copy )
!  print *, 'Integer_3D = ', Integer_3D
!  print *, 'Integer_3D_Copy = ', Integer_3D_Copy
!  print *
!
!  call Copy ( Integer_3D, [ 1, 1, 1 ], [ 1, 1, 1 ], [ 2, 2, 2 ], &
!              Integer_3D_Copy )
!  print *, 'Integer_3D = ', Integer_3D
!  print *, 'Integer_3D_Copy (Section) = ', Integer_3D_Copy
!  print *
!  print *

  !-- BigInteger

  BigInteger_1D = 1_KBI
!  call Clear ( BigInteger_1D_Copy )
  call Copy ( BigInteger_1D, BigInteger_1D_Copy )
  if ( all ( BigInteger_1D == BigInteger_1D_Copy ) ) then
    print *, 'Copy BigInteger_1D: Success'
  else
    print *, 'Copy BigInteger_1D: Fail'
    print *, 'BigInteger_1D = ', BigInteger_1D_Copy
    print *, 'BigInteger_1D_Copy = ', BigInteger_1D_Copy
  end if
  print *
  print *

!  BigInteger_2D = 1_KBI
!  call Clear ( BigInteger_2D_Copy )
!  print *, 'BigInteger_2D = ', BigInteger_2D
!  print *, 'BigInteger_2D_Copy = ', BigInteger_2D_Copy
!  print *

!  call Copy ( BigInteger_2D, BigInteger_2D_Copy )
!  print *, 'BigInteger_2D = ', BigInteger_2D
!  print *, 'BigInteger_2D_Copy = ', BigInteger_2D_Copy
!  print *
!  print *

!  BigInteger_3D = 1_KBI
!  call Clear ( BigInteger_3D_Copy )
!  print *, 'BigInteger_3D = ', BigInteger_3D
!  print *, 'BigInteger_3D_Copy = ', BigInteger_3D_Copy
!  print *

!  call Copy ( BigInteger_3D, BigInteger_3D_Copy )
!  print *, 'BigInteger_3D = ', BigInteger_3D
!  print *, 'BigInteger_3D_Copy = ', BigInteger_3D_Copy
!  print *
!  print *

  !-- BigInteger section

!  BigInteger_1D = 1_KBI
!  call Clear ( BigInteger_1D_Copy )
!  print *, 'BigInteger_1D = ', BigInteger_1D
!  print *, 'BigInteger_1D_Copy = ', BigInteger_1D_Copy
!  print *

!  call Copy ( BigInteger_1D, 1, 1, 2, BigInteger_1D_Copy )
!  print *, 'BigInteger_1D = ', BigInteger_1D
!  print *, 'BigInteger_1D_Copy (Section) = ', BigInteger_1D_Copy
!  print *
!  print *

!  BigInteger_2D = 1_KBI
!  call Clear ( BigInteger_2D_Copy )
!  print *, 'BigInteger_2D = ', BigInteger_2D
!  print *, 'BigInteger_2D_Copy = ', BigInteger_2D_Copy
!  print *

!  call Copy ( BigInteger_2D, [ 1, 1 ], [ 1, 1 ], [ 2, 2 ], BigInteger_2D_Copy )
!  print *, 'BigInteger_2D = ', BigInteger_2D
!  print *, 'BigInteger_2D_Copy (Section) = ', BigInteger_2D_Copy
!  print *
!  print *

!  BigInteger_3D = 1_KBI
!  call Clear ( BigInteger_3D_Copy )
!  print *, 'BigInteger_3D = ', BigInteger_3D
!  print *, 'BigInteger_3D_Copy = ', BigInteger_3D_Copy
!  print *

!  call Copy ( BigInteger_3D, [ 1, 1, 1 ], [ 1, 1, 1 ], [ 2, 2, 2 ], &
!              BigInteger_3D_Copy )
!  print *, 'BigInteger_3D = ', BigInteger_3D
!  print *, 'BigInteger_3D_Copy (Section) = ', BigInteger_3D_Copy
!  print *
!  print *

  !-- Real

  Real_1D = 1.0_KDR
  call Clear ( Real_1D_Copy )
  call Copy ( Real_1D, Real_1D_Copy )
  if ( all ( Real_1D == Real_1D_Copy ) ) then
    print *, 'Copy Real_1D: Success'
  else
    print *, 'Copy Real_1D: Fail'
    print *, 'Real_1D = ', Real_1D
    print *, 'Real_1D_Copy = ', Real_1D_Copy
  end if
  print *
  print *
  
  Real_1D = 1.0_KDR
  call AllocateDevice ( Real_1D, D_V )
  call AllocateDevice ( Real_1D_Copy, D_V_Copy )
  call AssociateHost  ( D_V,      Real_1D )
  call AssociateHost  ( D_V_Copy, Real_1D_Copy ) 
  
  call UpdateDevice ( Real_1D, D_V )
  call Clear ( Real_1D_Copy )
  call Clear ( Real_1D_Copy, UseDeviceOption = .true. )
  call Copy ( Real_1D, Real_1D_Copy, UseDeviceOption = .true. )
  call UpdateHost ( D_V_Copy, Real_1D_Copy )
  
  call DisassociateHost ( Real_1D_Copy )
  call DisassociateHost ( Real_1D )
  
  if ( all ( Real_1D == Real_1D_Copy ) ) then
    print *, 'Copy Real_1D UseDevice: Success'
  else
    print *, 'Copy Real_1D UseDevice: Fail'
    print *, 'Real_1D = ', Real_1D
    print *, 'Real_1D_Copy = ', Real_1D_Copy
  end if
  print *
  print *
  
  Real_2D = 1.0_KDR
  call Clear ( Real_2D_Copy )
  call Copy  ( Real_2D, Real_2D_Copy )
  if ( all ( Real_1D == Real_1D_Copy ) ) then
    print *, 'Copy Real_2D: Success'
  else
    print *, 'Copy Real_2D: Fail'
    print *, 'Real_2D = ', Real_1D
    print *, 'Real_2D_Copy = ', Real_1D_Copy
  end if
  print *
  print *
  
  Real_2D = 1.0_KDR
  call AllocateDevice ( Real_2D, D_V )
  call AllocateDevice ( Real_2D_Copy, D_V_Copy )
  call AssociateHost  ( D_V,      Real_2D )
  call AssociateHost  ( D_V_Copy, Real_2D_Copy ) 
  
  call UpdateDevice ( Real_2D, D_V )
  call Clear ( Real_2D_Copy )
  call Clear ( Real_2D_Copy, UseDeviceOption = .true. )
  call Copy ( Real_2D, Real_2D_Copy, UseDeviceOption = .true. )
  call UpdateHost ( D_V_Copy, Real_2D_Copy )
  
  call DisassociateHost ( Real_2D_Copy )
  call DisassociateHost ( Real_2D )
  
  if ( all ( Real_2D == Real_2D_Copy ) ) then
    print *, 'Copy Real_2D UseDevice: Success'
  else
    print *, 'Copy Real_2D UseDevice: Fail'
    print *, 'Real_2D = ', Real_2D
    print *, 'Real_2D_Copy = ', Real_2D_Copy
  end if
  print *
  print *
  
  Real_3D = 1.0_KDR
  call Clear ( Real_3D_Copy )
  call Copy ( Real_3D, Real_3D_Copy )
  if ( all ( Real_3D == Real_3D_Copy ) ) then
    print *, 'Copy Real_3D: Success'
  else
    print *, 'Copy Real_3D: Fail'
    print *, 'Real_3D = ', Real_3D
    print *, 'Real_3D_Copy = ', Real_3D_Copy
  end if
  print *
  print *
  
  Real_3D = 1.0_KDR
  call AllocateDevice ( Real_3D, D_V )
  call AllocateDevice ( Real_3D_Copy, D_V_Copy )
  call AssociateHost  ( D_V,      Real_3D )
  call AssociateHost  ( D_V_Copy, Real_3D_Copy ) 
  
  call UpdateDevice ( Real_3D, D_V )
  call Clear ( Real_3D_Copy )
  call Clear ( Real_3D_Copy, UseDeviceOption = .true. )
  call Copy ( Real_3D, Real_3D_Copy, UseDeviceOption = .true. )
  call UpdateHost ( D_V_Copy, Real_3D_Copy )
  
  call DisassociateHost ( Real_3D_Copy )
  call DisassociateHost ( Real_3D )
  
  if ( all ( Real_3D == Real_3D_Copy ) ) then
    print *, 'Copy Real_3D UseDevice: Success'
  else
    print *, 'Copy Real_3D UseDevice: Fail'
    print *, 'Real_3D = ', Real_3D
    print *, 'Real_3D_Copy = ', Real_3D_Copy
  end if
  print *
  print *
  
  return
  !-- Real section

!  Real_1D = 1.0_KDR
!  call Clear ( Real_1D_Copy )
!  print *, 'Real_1D = ', Real_1D
!  print *, 'Real_1D_Copy = ', Real_1D_Copy
!  print *

!  call Copy ( Real_1D, 1, 1, 2, Real_1D_Copy )
!  print *, 'Real_1D = ', Real_1D
!  print *, 'Real_1D_Copy (Section) = ', Real_1D_Copy
!  print *
!  print *

!  Real_2D = 1.0_KDR
!  call Clear ( Real_2D_Copy )
!  print *, 'Real_2D = ', Real_2D
!  print *, 'Real_2D_Copy = ', Real_2D_Copy
!  print *

!  call Copy ( Real_2D, [ 1, 1 ], [ 1, 1 ], [ 2, 2 ], Real_2D_Copy )
!  print *, 'Real_2D = ', Real_2D
!  print *, 'Real_2D_Copy (Section) = ', Real_2D_Copy
!  print *
!  print *

!  Real_3D = 1.0_KDR
!  call Clear ( Real_3D_Copy )
!  print *, 'Real_3D = ', Real_3D
!  print *, 'Real_3D_Copy = ', Real_3D_Copy
!  print *

!  call Copy ( Real_3D, [ 1, 1, 1 ], [ 1, 1, 1 ], [ 2, 2, 2 ], Real_3D_Copy )
!  print *, 'Real_3D = ', Real_3D
!  print *, 'Real_3D_Copy (Section) = ', Real_3D_Copy
!  print *
!  print *

  !-- Complex

!  Complex_1D = ( 1.0_KDR, 1.0_KDR )
!  call Clear ( Complex_1D_Copy )
!  print *, 'Complex_1D = ', Complex_1D
!  print *, 'Complex_1D_Copy = ', Complex_1D_Copy
!  print *

!  call Copy ( Complex_1D, Complex_1D_Copy )
!  print *, 'Complex_1D = ', Complex_1D
!  print *, 'Complex_1D_Copy = ', Complex_1D_Copy
!  print *
!  print *

!  Complex_2D = ( 1.0_KDR, 1.0_KDR )
!  call Clear ( Complex_2D_Copy )
!  print *, 'Complex_2D = ', Complex_2D
!  print *, 'Complex_2D_Copy = ', Complex_2D_Copy
!  print *

!  call Copy ( Complex_2D, Complex_2D_Copy )
!  print *, 'Complex_2D = ', Complex_2D
!  print *, 'Complex_2D_Copy = ', Complex_2D_Copy
!  print *
!  print *

!  Complex_3D = ( 1.0_KDR, 1.0_KDR )
!  call Clear ( Complex_3D_Copy )
!  print *, 'Complex_3D = ', Complex_3D
!  print *, 'Complex_3D_Copy = ', Complex_3D_Copy
!  print *

!  call Copy ( Complex_3D, Complex_3D_Copy )
!  print *, 'Complex_3D = ', Complex_3D
!  print *, 'Complex_3D_Copy = ', Complex_3D_Copy
!  print *
!  print *

  !-- Complex section

!  Complex_1D = ( 1.0_KDR, 1.0_KDR )
!  call Clear ( Complex_1D_Copy )
!  print *, 'Complex_1D = ', Complex_1D
!  print *, 'Complex_1D_Copy = ', Complex_1D_Copy
!  print *

!  call Copy ( Complex_1D, 1, 1, 2, Complex_1D_Copy )
!  print *, 'Complex_1D = ', Complex_1D
!  print *, 'Complex_1D_Copy (Section) = ', Complex_1D_Copy
!  print *
!  print *

!  Complex_2D = ( 1.0_KDR, 1.0_KDR )
!  call Clear ( Complex_2D_Copy )
!  print *, 'Complex_2D = ', Complex_2D
!  print *, 'Complex_2D_Copy = ', Complex_2D_Copy
!  print *

!  call Copy ( Complex_2D, [ 1, 1 ], [ 1, 1 ], [ 2, 2 ], Complex_2D_Copy )
!  print *, 'Complex_2D = ', Complex_2D
!  print *, 'Complex_2D_Copy (Section) = ', Complex_2D_Copy
!  print *
!  print *

!  Complex_3D = ( 1.0_KDR, 1.0_KDR )
!  call Clear ( Complex_3D_Copy )
!  print *, 'Complex_3D = ', Complex_3D
!  print *, 'Complex_3D_Copy = ', Complex_3D_Copy
!  print *

!  call Copy ( Complex_3D, [ 1, 1, 1 ], [ 1, 1, 1 ], [ 2, 2, 2 ], &
!              Complex_3D_Copy )
!  print *, 'Complex_3D = ', Complex_3D
!  print *, 'Complex_3D_Copy (Section) = ', Complex_3D_Copy
!  print *
!  print *

  !-- Character to String of C Char

!  Character_3 = 'Foo'
!  Character_1D = ''
!  print *, 'Character_3 = ', Character_3
!  print *, 'Character_1D = ', Character_1D
!  print *

!  call Copy ( Character_3, Character_1D )
!  print *, 'Character_3 = ', Character_3
!  print *, 'Character_1D ( 1 ) = ', Character_1D ( 1 )
!  print *, 'Character_1D ( 2 ) = ', Character_1D ( 2 )
!  print *, 'Character_1D ( 3 ) = ', Character_1D ( 3 )
!  print *, 'Character_1D ( 4 ) = ', Character_1D ( 4 )
!  print *, 'Character_1D = ', Character_1D
!  print *
!  print *
  
end program Copy_Command_Test
