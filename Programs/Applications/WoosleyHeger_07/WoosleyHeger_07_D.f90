program WoosleyHeger_07_D

  !-- WoosleyHeger_07_Deleptonization

  use GenASiS
  use WoosleyHeger_07_D__Form

  implicit none

  type ( WoosleyHeger_07_D_Form ), allocatable :: &
    WH

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07_D', DimensionalityOption = '1D' )

  allocate ( WH )
  call WH % Initialize ( PROGRAM_HEADER % Name )
  call WH % Evolve ( )
  deallocate ( WH )

  deallocate ( PROGRAM_HEADER )

end program WoosleyHeger_07_D
