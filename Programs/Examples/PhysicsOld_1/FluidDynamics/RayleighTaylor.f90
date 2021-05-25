program RayleighTaylor

  use GenASiS
  use RayleighTaylor_Form

  implicit none

  type ( RayleighTaylorForm ), allocatable :: &
    RT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RayleighTaylor', DimensionalityOption = '2D' )

  allocate ( RT )
  call RT % Initialize ( PROGRAM_HEADER % Name )
  call RT % Evolve ( )
  deallocate ( RT )

  deallocate ( PROGRAM_HEADER )

end program RayleighTaylor
