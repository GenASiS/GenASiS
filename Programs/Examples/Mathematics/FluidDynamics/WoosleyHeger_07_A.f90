program WoosleyHeger_07_A

  !-- WoosleyHeger_07_Adiabatic

  use Basics
  use WoosleyHeger_07_A__Form

  implicit none

  type ( WoosleyHeger_07_A_Form ), allocatable :: &
    FMT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'WoosleyHeger_07_A', DimensionalityOption = '1D' )

  allocate ( FMT )
  call FMT % Initialize ( PROGRAM_HEADER % Name )
  call FMT % Evolve ( )
  deallocate ( FMT )

  deallocate ( PROGRAM_HEADER )

end program WoosleyHeger_07_A
