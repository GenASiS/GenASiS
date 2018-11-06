program Geometry_ASC__Form_Test

  use Basics
  use Mathematics
  use Geometry_ASC__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Geometry_ASC__Form_Test'

  type ( Atlas_SC_Form ), allocatable :: &
    A_G, A_N, A_N_S
  type ( Geometry_ASC_Form ), allocatable :: &
    GA_G, GA_N, GA_N_S

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( A_G, A_N, A_N_S )
  call A_G % Initialize ( 'Atlas_SC_G', PROGRAM_HEADER % Communicator )
  call A_N % Initialize ( 'Atlas_SC_N', PROGRAM_HEADER % Communicator )
  call A_N_S % Initialize ( 'Atlas_SC_N_S', PROGRAM_HEADER % Communicator )
  call A_G % CreateChart ( )
  call A_N % CreateChart ( )
  call A_N_S % CreateChart ( )

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  allocate ( GA_G, GA_N, GA_N_S )
  call GA_G % Initialize &
         ( A_G, 'GALILEAN', NameShortOption = 'Geometry_G' )
  call GA_N % Initialize &
         ( A_N, 'NEWTONIAN', NameShortOption = 'Geometry_N', &
           GravitySolverTypeOption = 'MULTIPOLE' )
  call GA_N_S % Initialize &
         ( A_N_S, 'NEWTONIAN_STRESS', NameShortOption = 'Geometry_N_S', &
           GravitySolverTypeOption = 'MULTIPOLE' )
  deallocate ( GA_N_S, GA_N, GA_G )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( A_N_S, A_N, A_G )
  deallocate ( PROGRAM_HEADER )

end program Geometry_ASC__Form_Test
