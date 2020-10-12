module LaplacianMultipole_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use LaplacianMultipoleOld_1_ASC__Form
  use SetHomogeneousSphere_Command

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'LaplacianMultipole_ASC__Form_Test'
    
end module LaplacianMultipole_ASC__Form_Test__Form


program LaplacianMultipole_ASC__Form_Test

  use Basics
  use LaplacianMultipole_ASC__Form_Test__Form

  implicit none

!  type ( LaplacianMultipole_ASC__Form_Test_Form ), allocatable :: &
!    LMFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

!  allocate ( LMFT )
!  call LMFT % Initialize ( PROGRAM_HEADER % Name )
!  deallocate ( LMFT )

  deallocate ( PROGRAM_HEADER )

end program LaplacianMultipole_ASC__Form_Test
