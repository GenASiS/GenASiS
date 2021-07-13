program FluidBox_Form_Test

  use Basics
  use FluidBox_Form

  implicit none

  type ( FluidBoxForm ), allocatable :: &
    FB

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FluidBox_Form_Test', DimensionalityOption = '2D' )

  allocate ( FB )
  call FB % Initialize &
         ( FluidType = 'DUST', &
           GravitationType = 'GALILEO' )
  call FB % Show ( )
  deallocate ( FB )

  deallocate ( PROGRAM_HEADER )

end program FluidBox_Form_Test
