program RectangularExpansion_VLC

  !-- RectangularExpansion_VelocityLinearConstant

  use Basics
  use Mathematics
  use ImposedAdvection_VL__Form

  implicit none

  type ( ImposedAdvection_VL_Form ), allocatable :: &
    RE

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RectangularExpansion_VLC', DimensionalityOption = '1D' )

  allocate ( RE )
  call RE % Initialize &
         ( CoordinateSystem = 'RECTANGULAR', &
           AdvectionType = 'EXPANSION' )
  call RE % Evolve ( )
  deallocate ( RE )

  deallocate ( PROGRAM_HEADER )

end program RectangularExpansion_VLC
