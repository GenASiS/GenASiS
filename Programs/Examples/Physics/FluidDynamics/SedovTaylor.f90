program SedovTaylor

  use GenASiS
  use SedovTaylor_Form

  implicit none

  type ( SedovTaylorForm ), allocatable :: &
    ST

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SedovTaylor' )

  allocate ( ST )
  call ST % Initialize ( PROGRAM_HEADER % Name )
  call ST % Evolve ( )
  deallocate ( ST )

  deallocate ( PROGRAM_HEADER )

end program SedovTaylor
