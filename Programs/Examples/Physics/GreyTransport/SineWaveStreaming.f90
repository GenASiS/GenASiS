program SineWaveStreaming

  use Basics
  use SineWaveStreaming_Form

  implicit none

  type ( SineWaveStreamingForm ), allocatable :: &
    SWS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveStreaming', DimensionalityOption = '2D' )

  allocate ( SWS )
  call SWS % Initialize ( PROGRAM_HEADER % Name )
  call SWS % Evolve ( )
!  call SWS % ComputeError ( )
  deallocate ( SWS )

  deallocate ( PROGRAM_HEADER )

end program SineWaveStreaming
