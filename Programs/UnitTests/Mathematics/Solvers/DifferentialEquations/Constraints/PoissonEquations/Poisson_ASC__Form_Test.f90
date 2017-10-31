module Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Poisson_ASC__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Poisson_ASC__Form_Test'

end module Poisson_ASC__Form_Test__Form



program Poisson_ASC__Form_Test

  use Basics
  use Poisson_ASC__Form_Test__Form

  implicit none

!  type ( Poisson_ASC__Form_Test_Form ), allocatable :: &
!    LMFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

!  allocate ( LMFT )
!  call LMFT % Initialize ( PROGRAM_HEADER % Name )
!  deallocate ( LMFT )

  deallocate ( PROGRAM_HEADER )

end program Poisson_ASC__Form_Test
