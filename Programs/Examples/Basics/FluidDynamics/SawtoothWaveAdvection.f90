program SawtoothWaveAdvection

  use Basics
  use SawtoothWaveAdvection_Form

  implicit none

  type ( SawtoothWaveAdvectionForm ) :: &
    SWA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SawtoothWaveAdvection' )

  call SWA % Initialize ( )
  call SWA % Evolve ( )

  deallocate ( PROGRAM_HEADER )

end program SawtoothWaveAdvection
