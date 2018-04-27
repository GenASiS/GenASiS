program YahilLattimer

  use GenASiS
  use YahilLattimer_Form

  implicit none

  type ( YahilLattimerForm ), allocatable :: &
    YL

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'YahilLattimer', DimensionalityOption = '1D' )

  allocate ( YL )
  call YL % Initialize ( PROGRAM_HEADER % Name )
  call YL % Evolve ( )
!  call YL % ComputeError ( )
  deallocate ( YL )

  deallocate ( PROGRAM_HEADER )

end program YahilLattimer
