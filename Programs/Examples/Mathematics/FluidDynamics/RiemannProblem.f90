program RiemannProblem

  use Basics
  use RiemannProblem_Form

  implicit none

  type ( RiemannProblemForm ), allocatable :: &
    RP

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'RiemannProblem' )

  allocate ( RP )
  call RP % Initialize ( PROGRAM_HEADER % Name )
  call RP % Evolve ( )
  deallocate ( RP )

  deallocate ( PROGRAM_HEADER )

end program RiemannProblem
