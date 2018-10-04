program SpheroidCollapse

  use GenASiS
  use SpheroidCollapse_Form

  implicit none

  type ( SpheroidCollapseForm ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SpheroidCollapse', DimensionalityOption = '1D' )

  allocate ( SC )
  call SC % Initialize ( PROGRAM_HEADER % Name )
  call SC % Evolve ( )
!  call SC % ComputeError ( )
  deallocate ( SC )

  deallocate ( PROGRAM_HEADER )

end program SpheroidCollapse
