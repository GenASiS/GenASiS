program SineWaveStreaming_S

  use Basics
  use SineWaveStreaming_Form

  implicit none

  type ( SineWaveStreamingForm ), allocatable :: &
    SWS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveStreaming_S', DimensionalityOption = '2D_1D' )

  allocate ( SWS )
  call SWS % Initialize ( 'SPECTRAL', PROGRAM_HEADER % Name )
  call SWS % Evolve ( )
  call SWS % ComputeError ( )
  deallocate ( SWS )

  deallocate ( PROGRAM_HEADER )

end program SineWaveStreaming_S
