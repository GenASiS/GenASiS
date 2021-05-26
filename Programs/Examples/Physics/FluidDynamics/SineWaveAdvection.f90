program SineWaveAdvection

  use GenASiS
  use PlaneWave_Form
  use SineWave_Form

  implicit none

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveAdvection', DimensionalityOption = '2D' )

  allocate ( SineWaveForm :: PLANE_WAVE )
  associate ( SW  =>  PLANE_WAVE )

  call SW % Initialize ( )
  call SW % Evolve ( )
  call SW % ComputeError ( )

  end associate !-- SW
  deallocate ( PLANE_WAVE )

  deallocate ( PROGRAM_HEADER )

end program SineWaveAdvection
