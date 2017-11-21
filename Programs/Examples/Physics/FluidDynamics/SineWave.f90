program SineWave

  use GenASiS
  use SineWave_Form

  implicit none

  type ( SineWaveForm ), allocatable :: &
    SW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SineWave' )

  allocate ( SW )
  call SW % Initialize ( PROGRAM_HEADER % Name )
  call SW % Evolve ( )
  call SW % ComputeError ( )
  deallocate ( SW )

  deallocate ( PROGRAM_HEADER )

end program SineWave
