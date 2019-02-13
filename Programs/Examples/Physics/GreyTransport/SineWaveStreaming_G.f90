program SineWaveStreaming_G

  use Basics
  use SineWaveStreaming_G__Form

  implicit none

  type ( SineWaveStreaming_G_Form ), allocatable :: &
    SWS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveStreaming_G', DimensionalityOption = '2D' )

  allocate ( SWS )
  call SWS % Initialize ( PROGRAM_HEADER % Name )
  call SWS % Evolve ( )
  call SWS % ComputeError ( )
  deallocate ( SWS )

  deallocate ( PROGRAM_HEADER )

end program SineWaveStreaming_G
