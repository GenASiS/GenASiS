program MarshakWave

  use Basics
  use MarshakWave_Form

  implicit none

  type ( MarshakWaveForm ), allocatable :: &
    MW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'MarshakWave', DimensionalityOption = '1D_1D' )

  allocate ( MW )
  call MW % Initialize ( PROGRAM_HEADER % Name )
!  call MW % Evolve ( )
  deallocate ( MW )

  deallocate ( PROGRAM_HEADER )

end program MarshakWave
