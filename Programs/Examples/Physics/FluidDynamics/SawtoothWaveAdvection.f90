program SawtoothWaveAdvection

  use GenASiS
  use SawtoothWaveAdvection_Form

  implicit none

  type ( SawtoothWaveAdvectionForm ), allocatable :: &
    SWA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SawtoothWaveAdvection', DimensionalityOption = '2D' )

  allocate ( SWA )
  call SWA % Initialize ( PROGRAM_HEADER % Name )
  call SWA % Evolve ( )
  call SWA % ComputeError ( )
  deallocate ( SWA )

  deallocate ( PROGRAM_HEADER )

end program SawtoothWaveAdvection
