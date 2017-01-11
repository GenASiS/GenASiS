program SawtoothWave

  use Basics
  use SawtoothWave_Form

  implicit none

  type ( SawtoothWaveForm ), allocatable :: &
    SW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SawtoothWave' )

  allocate ( SW )
  call SW % Initialize ( PROGRAM_HEADER % Name )
  call SW % Evolve ( )
  call SW % ComputeError ( )
  deallocate ( SW )

  deallocate ( PROGRAM_HEADER )

end program SawtoothWave
