program WoosleyHeger_07

  use Basics
  use WoosleyHeger_07__Form

  implicit none

  type ( WoosleyHeger_07_Form ), allocatable :: &
    FMT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07', DimensionalityOption = '1D' )

  allocate ( FMT )
  call FMT % Initialize ( PROGRAM_HEADER % Name )
  call FMT % Evolve ( )
  deallocate ( FMT )

  deallocate ( PROGRAM_HEADER )

end program WoosleyHeger_07
