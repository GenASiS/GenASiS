program OppenheimerSnyder

  use GenASiS
  use OppenheimerSnyder_Form

  implicit none

  type ( OppenheimerSnyderForm ), allocatable :: &
    OS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'OppenheimerSnyder', DimensionalityOption = '1D' )

  allocate ( OS )
  call OS % Initialize ( PROGRAM_HEADER % Name )
  call OS % Evolve ( )
!  call OS % ComputeError ( )
  deallocate ( OS )

  deallocate ( PROGRAM_HEADER )

end program OppenheimerSnyder
