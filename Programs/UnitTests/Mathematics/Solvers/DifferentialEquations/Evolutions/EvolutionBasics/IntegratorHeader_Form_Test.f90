module IntegratorHeader_Form_Test__Form

  use Basics
  use IntegratorHeader_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'IntegratorHeader_Form_Test'

!   type, public :: IntegratorHeader_Form_Test_Form
!     type ( IntegratorHeaderForm ), allocatable :: &
!       IntegratorHeaderBase, &
!       IntegratorHeaderFiber
!   contains
!     procedure, public, pass :: &
!       Initialize
!     final :: &
!       Finalize
!   end type IntegratorHeader_Form_Test_Form

! contains

  
!   subroutine Initialize ( AHFT )

!     class ( IntegratorHeader_Form_Test_Form ), intent ( inout ), target :: &
!       AHFT

!     allocate ( AHFT % IntegratorHeaderBase )
!     associate ( A => AHFT % IntegratorHeaderBase )
!     call A % Initialize &
!            ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
!              iDimensionalityOption = 1 )
!     call A % Connectivity % Show ( )
!     call A % Show ( )
!     end associate !-- AH

!     allocate ( AHFT % IntegratorHeaderFiber )
!     associate ( A => AHFT % IntegratorHeaderFiber )
!     call A % Initialize &
!            ( 'Fiber', iDimensionalityOption = 2 )
!     call A % Connectivity % Show ( )
!     call A % Show ( )
!     end associate !-- AH

!   end subroutine Initialize


!   subroutine Finalize ( AHFT )

!     type ( IntegratorHeader_Form_Test_Form ), intent ( inout ) :: &
!       AHFT

!     deallocate ( AHFT % IntegratorHeaderFiber )
!     deallocate ( AHFT % IntegratorHeaderBase )

!   end subroutine Finalize


end module IntegratorHeader_Form_Test__Form



program IntegratorHeader_Form_Test

  use Basics
  use IntegratorHeader_Form_Test__Form
  
  implicit none

!  type ( IntegratorHeader_Form_Test_Form ), allocatable :: &
!    IHFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
!  allocate ( IHFT )
!  call IHFT % Initialize ( )
!  deallocate ( IHFT )

  deallocate ( PROGRAM_HEADER )

end program IntegratorHeader_Form_Test
