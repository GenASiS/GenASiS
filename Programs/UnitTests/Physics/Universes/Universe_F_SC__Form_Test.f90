program Universe_F_SC__Form_Test

  !-- Universe_Fluid_SymmetricCurvilinear__Form_Test

  use Basics
  use Universe_F_SC__Form

  implicit none

  type ( Universe_F_SC_Form ), allocatable :: &
    U

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Universe_F_SC__Form_Test', DimensionalityOption = '2D' )

  allocate ( U )
  call U % Initialize &
         ( FluidType = 'IDEAL', &
           Name = 'Universe', &
           RadiusMax = 10.0_KDR )
  call U % Show ( )
  deallocate ( U )

  deallocate ( PROGRAM_HEADER )

end program Universe_F_SC__Form_Test
