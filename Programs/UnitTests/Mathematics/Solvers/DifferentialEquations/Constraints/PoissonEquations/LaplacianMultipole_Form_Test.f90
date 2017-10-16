module LaplacianMultipole_Form_Test__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'LaplacianMultipole_Form_Test'
    
end module LaplacianMultipole_Form_Test__Form


program LaplacianMultipole_Form_Test

  use Basics
  use LaplacianMultipole_Form_Test__Form

  implicit none

!  type ( LaplacianMultipole_Form_Test_Form ), allocatable :: &
!    LMFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

!  allocate ( LMFT )
!  call LMFT % Initialize ( PROGRAM_HEADER % Name )
!  deallocate ( LMFT )

  deallocate ( PROGRAM_HEADER )

end program LaplacianMultipole_Form_Test

