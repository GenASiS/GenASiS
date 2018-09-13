program SASI

  use GenASiS
  use SASI_Form

  implicit none

  type ( SASIForm ), allocatable :: &
    SF

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SASI', DimensionalityOption = '2D' )

  allocate ( SF )
  call SF % Initialize ( PROGRAM_HEADER % Name )
  call SF % Evolve ( )
  deallocate ( SF )

  deallocate ( PROGRAM_HEADER )

end program SASI
