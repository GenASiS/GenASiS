program FluidCentralExcision_Form_Test

  use Basics
  use FluidCentralExcision_Form

  implicit none

  type ( FluidCentralExcisionForm ), allocatable :: &
    FCE

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FluidCentralExcision_Form_Test', DimensionalityOption = '2D' )

  allocate ( FCE )
  call FCE % Initialize &
         ( PROGRAM_HEADER % Name, &
           FluidType = 'IDEAL', &
           GeometryType = 'NEWTONIAN' )
  deallocate ( FCE )

  deallocate ( PROGRAM_HEADER )

end program FluidCentralExcision_Form_Test
