program FluidBox_Form_Test

  use Basics
  use FluidBox_Form

  implicit none

  type ( FluidBoxForm ), allocatable :: &
    FB

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'FluidBox_Form_Test' )

  allocate ( FB )
  call FB % Initialize &
         ( PROGRAM_HEADER % Name, &
           FluidType = 'DUST', &
           GeometryType = 'GALILEAN' )
  deallocate ( FB )

  deallocate ( PROGRAM_HEADER )

end program FluidBox_Form_Test
