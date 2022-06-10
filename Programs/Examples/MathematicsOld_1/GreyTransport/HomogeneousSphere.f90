program HomogeneousSphere

  use Basics
  use HomogeneousSphere_Form

  implicit none

  type ( HomogeneousSphereForm ), allocatable :: &
    HS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'HomogeneousSphere' )

  allocate ( HS )
  call HS % Initialize ( PROGRAM_HEADER % Name )
  call HS % Evolve ( )
  call HS % ComputeError ( ) 
  deallocate ( HS )

  deallocate ( PROGRAM_HEADER )

end program HomogeneousSphere
