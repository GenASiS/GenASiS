program HomogeneousSpheroid

  use Basics
  use HomogeneousSpheroid_Form
  
  implicit none

  type ( HomogeneousSpheroidForm ), allocatable :: &
    HS
  integer ( KDI ) :: &
    N_EQUATIONS = 3
  character ( LDL ), dimension ( : ), allocatable :: &
    VARIABLE
  real ( KDR ) :: &
    Density
  real ( KDR ), dimension ( : ), allocatable :: &
    SemiMajor, &
    Eccentricity

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'HomogeneousSphereoid_Test' )

  allocate &
    ( VARIABLE    ( N_EQUATIONS ), &
      SemiMajor   ( N_EQUATIONS ), &
      Eccentricity ( N_EQUATIONS ) )

  VARIABLE = [ 'OblateHomogeneousSpheroid_1 ', &
               'OblateHomogeneousSpheroid_2 ', &
               'ProlateHomogeneousSpheroid_1' ]

  allocate ( HS )

  call HS % Initialize ( PROGRAM_HEADER % Name, N_EQUATIONS, VARIABLE )

  SemiMajor = HS % Atlas % Chart % MaxCoordinate ( 1 ) / 2
  call PROGRAM_HEADER % GetParameter ( SemiMajor, 'SemiMajor' )

  Eccentricity = sqrt ( 1.0_KDR &
                         - [ 0.7_KDR ** 2, 0.2_KDR ** 2, 0.2_KDR ** 2 ] )
  call PROGRAM_HEADER % GetParameter ( Eccentricity, 'Eccentricity' )

  Density = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  )
  call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

  call HS % SetHomogeneousSpheroid ( Density, SemiMajor, Eccentricity )

  call HS % Compute ( )
  
  deallocate ( HS )

  deallocate ( PROGRAM_HEADER )

end program HomogeneousSpheroid