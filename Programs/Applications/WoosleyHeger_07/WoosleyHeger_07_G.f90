program WoosleyHeger_07_G

  !-- WoosleyHeger_07_Grey

  use GenASiS
  use WoosleyHeger_07_RM__Form

  implicit none

  type ( WoosleyHeger_07_RM_Form ), allocatable :: &
    WH

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07_G', DimensionalityOption = '1D' )

  allocate ( WH )
  call WH % Initialize ( 'GREY', PROGRAM_HEADER % Name )
  call WH % Evolve ( )
  deallocate ( WH )

  deallocate ( PROGRAM_HEADER )

end program WoosleyHeger_07_G
