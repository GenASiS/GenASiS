program Universe_F_B__Form_Test

  !-- Universe_Fluid_Box__Form_Test

  use Basics
  use Universe_F_B__Form

  implicit none

  type ( Universe_F_B_Form ), allocatable :: &
    U

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Universe_F_B__Form_Test', DimensionalityOption = '2D' )

  allocate ( U )
  call U % Initialize &
         ( FluidType = 'DUST', &
           GravitationType = 'GALILEO' )
  call U % Show ( )
  deallocate ( U )

  deallocate ( PROGRAM_HEADER )

end program Universe_F_B__Form_Test
