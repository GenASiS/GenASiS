program CondensedSphere

  use Basics
  use CondensedSphere_Form
  
  implicit none

  type ( CondensedSphereForm ), allocatable :: &
    CS
  integer ( KDI ) :: &
    N_EQUATIONS = 3
  character ( LDL ), dimension ( : ), allocatable :: &
    VARIABLE
  real ( KDR ) :: &
    Density
  real ( KDR ), dimension ( : ), allocatable :: &
    SR, &
    rr

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'CondensedSphere_Test' )

  allocate &
    ( VARIABLE ( N_EQUATIONS ), &
      SR       ( N_EQUATIONS ), &
      rr       ( N_EQUATIONS ) )

  VARIABLE = [ 'CondensedSphere_1           ', &
               'CondensedSphere_2           ', &
               'CondensedSphere_3           ' ]

  allocate ( CS )

  call CS % Initialize ( PROGRAM_HEADER % Name, N_EQUATIONS, VARIABLE )

  SR = [ 0.1_KDR, 0.5_KDR, 0.5_KDR ] * CS % Atlas % Chart % MaxCoordinate ( 1 )
  call PROGRAM_HEADER % GetParameter ( SR, 'SR' )

  rr = [ 0.05_KDR, 0.1_KDR, 0.3_KDR ] * CS % Atlas % Chart % MaxCoordinate ( 1 )
  call PROGRAM_HEADER % GetParameter ( rr, 'rr' )

  Density = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  )
  call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

  call CS % SetCondensedSphere ( Density, SR, rr )

  call CS % Compute ( )
  
  deallocate ( CS )

  deallocate ( PROGRAM_HEADER )

end program CondensedSphere
