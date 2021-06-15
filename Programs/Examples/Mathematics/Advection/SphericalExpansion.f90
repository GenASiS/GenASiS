program SphericalExpansion

  use Basics
  use Mathematics
  use LinearAdvection_Form

  implicit none

  type ( LinearAdvectionForm ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SphericalExpansion', DimensionalityOption = '1D' )

  allocate ( SC )
  call SC % Initialize &
         ( CoordinateSystem = 'SPHERICAL', &
           AdvectionType = 'EXPANSION' )
  call SC % Evolve ( )
  deallocate ( SC )

  deallocate ( PROGRAM_HEADER )

end program SphericalExpansion
