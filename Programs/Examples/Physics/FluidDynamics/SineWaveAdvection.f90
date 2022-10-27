program SineWaveAdvection

  use GenASiS
  use SineWaveAdvection_Form

  implicit none

  type ( SineWaveAdvectionForm ), allocatable :: &
    SWA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveAdvection', DimensionalityOption = '2D' )

  allocate ( SWA )
  call SWA % Initialize ( PROGRAM_HEADER % Name )
  call SWA % Evolve ( )
  call SWA % ComputeError ( )
  deallocate ( SWA )

  deallocate ( PROGRAM_HEADER )

end program SineWaveAdvection
