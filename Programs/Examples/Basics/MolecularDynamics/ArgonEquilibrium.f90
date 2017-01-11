program ArgonEquilibrium

  use Basics
  use LennardJonesDynamics_Form

  implicit none

  type ( LennardJonesDynamicsForm ) :: &
    LJD

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'ArgonEquilibrium' )

  call LJD % Initialize ( )
  call LJD % Evolve ( )

  deallocate ( PROGRAM_HEADER )

end program ArgonEquilibrium
