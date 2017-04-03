program WoosleyHeger_07_G

  !-- WoosleyHeger_07_Grey

  use Basics
  use WoosleyHeger_07_G__Form

  implicit none

  type ( WoosleyHeger_07_G_Form ), allocatable :: &
    FMT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07_G', DimensionalityOption = '1D' )

  allocate ( FMT )
  call FMT % Initialize ( PROGRAM_HEADER % Name )
!  call FMT % Evolve ( )
  deallocate ( FMT )

  deallocate ( PROGRAM_HEADER )

end program WoosleyHeger_07_G
