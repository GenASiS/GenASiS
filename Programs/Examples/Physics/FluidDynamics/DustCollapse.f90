program DustCollapse

  use GenASiS
  use DustCollapse_Form

  implicit none

  type ( DustCollapseForm ), allocatable :: &
    DC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'DustCollapse' )

  allocate ( DC )
  call DC % Initialize ( PROGRAM_HEADER % Name )
  call DC % Evolve ( )
!  call DC % ComputeError ( )
  deallocate ( DC )

  deallocate ( PROGRAM_HEADER )

end program DustCollapse
