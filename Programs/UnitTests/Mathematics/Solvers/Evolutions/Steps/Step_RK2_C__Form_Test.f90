program Step_RK2_C__Form_Test

  use Basics
  use DensityWaveStep_Form

  implicit none

  type ( DensityWaveStepForm ), allocatable :: &
    DW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'Step_RK2_C__Form_Test' )
    
  allocate ( DW )
  call DW % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( DW )

  deallocate ( PROGRAM_HEADER )

end program Step_RK2_C__Form_Test
