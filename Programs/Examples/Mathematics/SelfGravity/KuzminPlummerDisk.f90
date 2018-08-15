program KuzminPlummerDisk

  use Basics
  use KuzminPlummerDisk_Form
  
  implicit none

  type ( KuzminPlummerForm ), allocatable :: &
    KP
  integer ( KDI ) :: &
    N_EQUATIONS = 4
  character ( LDL ), dimension ( : ), allocatable :: &
    VARIABLE
  real ( KDR ) :: &
    TotalMass
  real ( KDR ), dimension ( : ), allocatable :: &
    a, &
    ratio

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'KuzminPlummerDisk_Test' )

  allocate &
    ( VARIABLE ( N_EQUATIONS ), &
      a        ( N_EQUATIONS ), &
      ratio    ( N_EQUATIONS ) )

  VARIABLE = [ 'KuzminPlummerDisk_1  ', &
               'KuzminPlummerDisk_2  ', &
               'KuzminPlummerDisk_3  ', &
               'KuzminPlummerDisk_4  ' ]

  allocate ( KP )

  call KP % Initialize ( PROGRAM_HEADER % Name, N_EQUATIONS, VARIABLE )

  a = 1.0_KDR
  call PROGRAM_HEADER % GetParameter ( a, 'a' )

  ratio = [ 0.2_KDR, 1.0_KDR, 5.0_KDR, 0.25_KDR ]
  call PROGRAM_HEADER % GetParameter ( ratio, 'ratio' )
  
  TotalMass = 1.0_KDR 
  call PROGRAM_HEADER % GetParameter ( TotalMass, 'TotalMass' )

  call KP % SetKuzminPlummer ( TotalMass, a, ratio )

  call KP % Compute ( NormalizeSolutionOption = .true. )
  
  deallocate ( KP )

  deallocate ( PROGRAM_HEADER )

end program KuzminPlummerDisk
