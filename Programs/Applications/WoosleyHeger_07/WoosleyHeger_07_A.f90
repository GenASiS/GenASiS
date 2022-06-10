program WoosleyHeger_07_A

  !-- WoosleyHeger_07_Adiabatic

  use GenASiS
  use WoosleyHeger_07_A__Form

  implicit none

  type ( WoosleyHeger_07_A_Form ), allocatable :: &
    WH

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07_A', DimensionalityOption = '1D' )

  allocate ( WH )
  call WH % Initialize ( PROGRAM_HEADER % Name )
  call WH % Evolve ( )
  deallocate ( WH )

  deallocate ( PROGRAM_HEADER )

end program WoosleyHeger_07_A
