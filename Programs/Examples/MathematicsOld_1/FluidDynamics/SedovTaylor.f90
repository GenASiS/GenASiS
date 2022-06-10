program SedovTaylor

  use Basics
  use SedovTaylor_Form

  implicit none

  type ( SedovTaylorForm ), allocatable :: &
    ST

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SedovTaylor', DimensionalityOption = '2D' )

  allocate ( ST )
  call ST % Initialize ( PROGRAM_HEADER % Name )
  call ST % Evolve ( )
  deallocate ( ST )

  deallocate ( PROGRAM_HEADER )

end program SedovTaylor
