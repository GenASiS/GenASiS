program SineWaveAdvection

  use Basics
  use SineWaveAdvection_Form

  implicit none

  type ( SineWaveAdvectionForm ) :: &
    SWA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SineWaveAdvection' )

  call SWA % Initialize ( )
  call SWA % Evolve ( )

  deallocate ( PROGRAM_HEADER )

end program SineWaveAdvection
