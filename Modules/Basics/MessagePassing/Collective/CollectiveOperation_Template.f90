!-- CollectiveOperationTemplate provides an abstraction commonly shared by 
!   higher-level, concrete type object with specified datatype for handling
!   collective operations.

module CollectiveOperation_Template

  use MPI
  use Specifiers
  use MessagePassingBasics

  implicit none
  private

  integer ( KDI ), public, parameter :: &
    UNSET = - huge ( 1_KDI )

  type, public, abstract :: CollectiveOperationTemplate
    integer ( KDI ) :: &
      Root = UNSET, &
      Error
    integer ( KDI ), dimension ( : ), allocatable :: &
      nIncoming, &
      nOutgoing
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
  end type CollectiveOperationTemplate
  
contains

  
  subroutine InitializeTemplate ( CO, C, nOutgoing, nIncoming, RootOption )
  
    class ( CollectiveOperationTemplate ), intent ( inout ) :: &
      CO
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ), dimension ( : ) :: &
      nOutgoing, &
      nIncoming
    integer ( KDI ), intent ( in ), optional :: &
      RootOption
  
    CO % Communicator => C

    if ( present ( RootOption ) ) CO % Root = RootOption

!-- FIXME: NAG 5.3.1 complains about this sourced allocation, even though
!          this feature is supposed to be supported according to the manual
!    allocate ( CO % nIncoming, source = nIncoming )
!    allocate ( CO % nOutgoing, source = nOutgoing )
    allocate ( CO % nIncoming ( size ( nIncoming ) ) )
    allocate ( CO % nOutgoing ( size ( nOutgoing ) ) )  
    CO % nIncoming = nIncoming
    CO % nOutgoing = nOutgoing

  end subroutine InitializeTemplate


end module CollectiveOperation_Template
