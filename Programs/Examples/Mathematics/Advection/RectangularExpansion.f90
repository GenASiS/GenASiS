program RectangularExpansion

  use Basics
  use Mathematics
  use LinearAdvection_Form

  implicit none

  type ( LinearAdvectionForm ), allocatable :: &
    RC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RectangularExpansion', DimensionalityOption = '1D' )

  allocate ( RC )
  call RC % Initialize &
         ( CoordinateSystem = 'RECTANGULAR', &
           AdvectionType = 'EXPANSION' )
  call RC % Evolve ( )
  deallocate ( RC )

  deallocate ( PROGRAM_HEADER )

end program RectangularExpansion
