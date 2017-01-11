program DensityWaveIntegrator_Form_Test

  use Basics
  use DensityWaveIntegrator_Form

  implicit none

  type ( DensityWaveIntegratorForm ), allocatable :: &
    DW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'DensityWaveIntegrator_Form_Test' )
    
  allocate ( DW )
  call DW % Initialize ( PROGRAM_HEADER % Name )
  call DW % Evolve ( )
  call DW % ComputeError ( )
  deallocate ( DW )

  deallocate ( PROGRAM_HEADER )

end program DensityWaveIntegrator_Form_Test
