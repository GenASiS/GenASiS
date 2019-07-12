program RiemannProblem

  use Basics
  use RiemannProblem_Form

  implicit none

  type ( RiemannProblemForm ) :: &
    RP

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'RiemannProblem' )

  call RP % Initialize ( )
  call RP % Evolve ( )
  
  deallocate ( PROGRAM_HEADER )

end program RiemannProblem
