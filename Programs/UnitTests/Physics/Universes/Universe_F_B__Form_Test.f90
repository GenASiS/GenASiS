program Universe_F_B__Form_Test

  !-- Universe_Fluid_Box__Form_Test

  use Basics
  use Universe_F_B__Form

  implicit none

  type ( Universe_F_B_Form ), allocatable :: &
    FB

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Universe_F_B__Form_Test', DimensionalityOption = '2D' )

  allocate ( FB )
  call FB % Initialize &
         ( FluidType = 'DUST', &
           GravitationType = 'GALILEO' )
  call FB % Show ( )
  deallocate ( FB )

  deallocate ( PROGRAM_HEADER )

end program Universe_F_B__Form_Test
