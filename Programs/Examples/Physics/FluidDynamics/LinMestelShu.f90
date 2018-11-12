program LinMestelShu

  use GenASiS
  use LinMestelShu_Form

  implicit none

  type ( LinMestelShuForm ), allocatable :: &
    LMS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'LinMestelShu', DimensionalityOption = '2D' )

  allocate ( LMS )
  call LMS % Initialize ( PROGRAM_HEADER % Name )
  call LMS % Evolve ( )
  call LMS % ComputeError ( )
  deallocate ( LMS )

  deallocate ( PROGRAM_HEADER )

end program LinMestelShu
