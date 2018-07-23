program TwoPowerDensity

  use Basics
  use TwoPowerDensity_Form
  
  implicit none

  type ( TwoPowerDensityForm ), allocatable :: &
    TP
  integer ( KDI ) :: &
    N_EQUATIONS = 3
  character ( LDL ), dimension ( : ), allocatable :: &
    VARIABLE
  real ( KDR ) :: &
    Density
  real ( KDR ), dimension ( : ), allocatable :: &
    a, &
    alpha, &
    beta

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'TwoPowerDensity_Test' )

  allocate &
    ( VARIABLE ( N_EQUATIONS ), &
      alpha    ( N_EQUATIONS ), &
      beta     ( N_EQUATIONS ), &
      a        ( N_EQUATIONS ) )

  VARIABLE = [ 'JaffeModel                  ', &
               'HernquistModel              ', &
               'NFW_Model                   ' ]

  allocate ( TP )

  call TP % Initialize ( PROGRAM_HEADER % Name, N_EQUATIONS, VARIABLE )

  alpha = [ 2.0_KDR, 1.0_KDR, 1.0_KDR ] 
  call PROGRAM_HEADER % GetParameter ( alpha, 'alpha' )

  beta  = [ 4.0_KDR, 4.0_KDR, 3.0_KDR ] 
  call PROGRAM_HEADER % GetParameter ( beta, 'beta' )

  a = 0.1_KDR
  call PROGRAM_HEADER % GetParameter ( a, 'a' ) 

  Density = 1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  )
  call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

  call TP % SetTwoPowerDensity ( Density, a, alpha, beta )

  call TP % Compute ( )
  
  deallocate ( TP )

  deallocate ( PROGRAM_HEADER )

end program TwoPowerDensity
