program MarshakWave_S

  use Basics
  use MarshakWave_Form

  implicit none

  type ( MarshakWaveForm ), allocatable :: &
    MW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'MarshakWave_S', DimensionalityOption = '1D' )

  allocate ( MW )
  call MW % Initialize ( 'SPECTRAL', PROGRAM_HEADER % Name )
  call MW % Evolve ( )
  deallocate ( MW )

  deallocate ( PROGRAM_HEADER )

end program MarshakWave_S
