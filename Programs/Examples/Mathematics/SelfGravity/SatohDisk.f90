program SatohDisk

  use Basics
  use SatohDisk_Form
  
  implicit none

  type ( SatohForm ), allocatable :: &
    S
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
  call PROGRAM_HEADER % Initialize ( 'SatohDisk_Test' )

  allocate &
    ( VARIABLE ( N_EQUATIONS ), &
      a        ( N_EQUATIONS ), &
      ratio    ( N_EQUATIONS ) )

  VARIABLE = [ 'SatohDisk_1  ', &
               'SatohDisk_2  ', &
               'SatohDisk_3  ', &
               'SatohDisk_4  ' ]

  allocate ( S )

  call S % Initialize ( PROGRAM_HEADER % Name, N_EQUATIONS, VARIABLE )

  a = 1.0_KDR
  call PROGRAM_HEADER % GetParameter ( a, 'a' )

  ratio = [ 0.2_KDR, 1.0_KDR, 5.0_KDR, 0.25_KDR ]
  call PROGRAM_HEADER % GetParameter ( ratio, 'ratio' )
  
  TotalMass = 1.0_KDR 
  call PROGRAM_HEADER % GetParameter ( TotalMass, 'TotalMass' )

  call S % SetSatoh ( TotalMass, a, ratio )

  call S % Compute ( NormalizeSolutionOption = .true. )
  
  deallocate ( S )

  deallocate ( PROGRAM_HEADER )

end program SatohDisk
