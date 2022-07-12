program Universe_F_CE__Form_Test

  !-- Universe_Fluid_CentralExcision__Form_Test

  use Basics
  use Mathematics
  use Fluids
  use Universe_F_CE__Form

  implicit none

  type ( TimerForm ), pointer :: &
    T_A, &
    T_W
  type ( Universe_F_CE_Form ), allocatable, target :: &
    U

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Universe_F_CE__Form_Test', DimensionalityOption = '2D' )

  allocate ( U )
  call U % Initialize &
         ( FluidType = 'DUST', &
           GravitationType = 'NEWTON_CM', &
           DimensionlessOption = .true., &
           CentralMassOption = 1.0_KDR )
  call U % Show ( )

  deallocate ( U )
  deallocate ( PROGRAM_HEADER )

end program Universe_F_CE__Form_Test
