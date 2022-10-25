program SawtoothWaveAdvection

  use GenASiS
  use SawtoothWave_Form

  implicit none

  type ( SawtoothWaveForm ), allocatable :: &
    SW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SawtoothWaveAdvection', DimensionalityOption = '2D' )

  allocate ( SW )
  call SW % Initialize ( Name = PROGRAM_HEADER % Name )
  call SW % Evolve ( )
  call SW % ComputeError ( )
  deallocate ( SW )

  deallocate ( PROGRAM_HEADER )

end program SawtoothWaveAdvection
