program LogarithmicPotential

  use Basics
  use LogarithmicPotential_Form
  
  implicit none

  type ( LogarithmicForm ), allocatable :: &
    L
  integer ( KDI ) :: &
    N_EQUATIONS = 3
  character ( LDL ), dimension ( : ), allocatable :: &
    VARIABLE
  real ( KDR ), dimension ( : ), allocatable :: &
    v0, &
    rho_c, &
    q_phi

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'LogarithmicPotential_Test' )

  allocate &
    ( VARIABLE ( N_EQUATIONS ), &
      v0       ( N_EQUATIONS ), &
      rho_c    ( N_EQUATIONS ), &
      q_phi    ( N_EQUATIONS ) )

  VARIABLE = [ 'LogarithmicPotential_1  ', &
               'LogarithmicPotential_2  ', &
               'LogarithmicPotential_3  ' ]

  allocate ( L )

  call L % Initialize ( PROGRAM_HEADER % Name, N_EQUATIONS, &
                        VARIABLE, RadiusMaxOption = 10.0_KDR )

  v0 = 1.0_KDR
  call PROGRAM_HEADER % GetParameter ( v0, 'v0' )

  rho_c = 0.1_KDR 
  call PROGRAM_HEADER % GetParameter ( rho_c, 'rho_c' )

  q_phi = [ 1.0_KDR/ sqrt ( 2.0_KDR ), 0.75_KDR, 0.8_KDR ]
  call PROGRAM_HEADER % GetParameter ( q_phi, 'q_phi' )

  call L % SetLogarithmic ( v0, rho_c, q_phi )

  call L % Compute &
         ( ShiftSolutionOption = .true. )
  
  deallocate ( L )

  deallocate ( PROGRAM_HEADER )

end program LogarithmicPotential
