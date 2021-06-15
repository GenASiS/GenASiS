program SphericalContraction

  use Basics
  use Mathematics
  use LinearAdvection_Form

  implicit none

  type ( LinearAdvectionForm ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SphericalContraction', DimensionalityOption = '1D' )

  allocate ( SC )
  call SC % Initialize &
         ( CoordinateSystem = 'SPHERICAL', &
           AdvectionType = 'CONTRACTION' )
  call SC % Evolve ( )
  deallocate ( SC )

  deallocate ( PROGRAM_HEADER )

end program SphericalContraction
