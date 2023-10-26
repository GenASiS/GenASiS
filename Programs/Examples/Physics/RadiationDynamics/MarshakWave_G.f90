program MarshakWave_G

  use Basics
  use MarshakWave_Form

  implicit none

  type ( MarshakWaveForm ), allocatable :: &
    MW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'MarshakWave_G', DimensionalityOption = '1D' )

  allocate ( MW )
  call MW % Initialize ( 'GREY', PROGRAM_HEADER % Name )
  call MW % Evolve ( )
  deallocate ( MW )

  deallocate ( PROGRAM_HEADER )

end program MarshakWave_G
