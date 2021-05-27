program SineWaveAdvection

  use GenASiS
  use PlaneWave_Form
  use SineWave_Form

  implicit none

  type ( SineWaveForm ), allocatable :: &
    SW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveAdvection', DimensionalityOption = '2D' )

  allocate ( SW )
  call SW % Initialize ( )
  call SW % Evolve ( )
  call SW % ComputeError ( )
  deallocate ( SW )

  deallocate ( PROGRAM_HEADER )

end program SineWaveAdvection
