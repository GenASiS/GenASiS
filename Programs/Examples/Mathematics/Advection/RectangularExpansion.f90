program RectangularExpansion

  use Basics
  use Mathematics
  use LinearAdvection_Form

  implicit none

  type ( LinearAdvectionForm ), allocatable :: &
    RE

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RectangularExpansion', DimensionalityOption = '1D' )

  allocate ( RE )
  call RE % Initialize &
         ( CoordinateSystem = 'RECTANGULAR', &
           AdvectionType = 'EXPANSION' )
  call RE % Evolve ( )
  deallocate ( RE )

  deallocate ( PROGRAM_HEADER )

end program RectangularExpansion
