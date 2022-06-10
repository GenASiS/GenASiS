program RectangularContraction_VLC

  !-- RectangularContraction_VelocityLinearConstant

  use Basics
  use Mathematics
  use ImposedAdvection_VL__Form

  implicit none

  type ( ImposedAdvection_VL_Form ), allocatable :: &
    RC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RectangularContraction_VLC', DimensionalityOption = '1D' )

  allocate ( RC )
  call RC % Initialize &
         ( CoordinateSystem = 'RECTANGULAR', &
           AdvectionType = 'CONTRACTION' )
  call RC % Evolve ( )
  deallocate ( RC )

  deallocate ( PROGRAM_HEADER )

end program RectangularContraction_VLC
