program SineWave

  use Basics
  use SineWave_Form

  implicit none

  type ( SineWaveForm ), allocatable :: &
    SW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SineWave', DimensionalityOption = '2D' )

  allocate ( SW )
  call SW % Initialize ( PROGRAM_HEADER % Name )
  call SW % Evolve ( )
  call SW % ComputeError ( )
  deallocate ( SW )

  deallocate ( PROGRAM_HEADER )

end program SineWave
