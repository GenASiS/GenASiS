program ClusterFormation
  
  use Basics
  use GravitationalDynamics_Form

  implicit none

  type ( GravitationalDynamicsForm ) :: &
    GD

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'ClusterFormation' )

  call GD % Initialize ( )
  call GD % Evolve ( )

  deallocate ( PROGRAM_HEADER )

end program ClusterFormation
