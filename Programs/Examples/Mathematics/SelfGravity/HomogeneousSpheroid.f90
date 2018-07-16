program HomogeneousSpheroid

  use Basics
  use HomogeneousSpheroid_Form
  
  implicit none

  type ( HomogeneousSpheroidForm ), allocatable :: &
    HS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'HomogeneousSphereoid_Test' )

  allocate ( HS )
  call HS % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( HS )

  deallocate ( PROGRAM_HEADER )

end program HomogeneousSpheroid
