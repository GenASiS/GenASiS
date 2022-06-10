program SphericalExpansion

  !-- SphericalContraction_VelocityLinearConstant

  use Basics
  use Mathematics
  use ImposedAdvection_VL__Form

  implicit none

  type ( ImposedAdvection_VL_Form ), allocatable :: &
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
