!-- MessageIncoming_1D_C_Form provides the concrete extension of 
!   Message_1D_Template for complex datatype to handle receiving array of 
!   messages.

module MessageIncoming_1D_C__Form

  use MPI
  use VariableManagement
  use MessagePassingBasics
  use Message_Template
  use MessageIncoming_C__Form
  use Message_1D__Template 
 
  implicit none
  private

  type, public, extends ( Message_1D_Template ) :: MessageIncoming_1D_C_Form
    type ( MessageIncoming_C_Form ), dimension ( : ), allocatable :: &
      Message
  contains
    procedure, public, pass :: &
      Initialize => InitializeAllocate
    procedure, public, pass :: &
      Receive
    final :: &
      Finalize
  end type MessageIncoming_1D_C_Form

contains


  subroutine InitializeAllocate ( M_1D, C, Tag, Rank, nElements )

    class ( MessageIncoming_1D_C_Form ), intent ( inout ), target :: &
      M_1D
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      Tag, &
      Rank, &
      nElements

    integer ( KDI ) :: &
      nMessages
    integer ( KDI ) :: &
      iM    !-- iMessage
    type ( MessageIncoming_C_Form ), pointer :: &
      M      

    nMessages = size ( Tag )

    allocate ( M_1D % Message ( nMessages ) )

    M_1D % nMessages = size ( Tag )
    M_1D % MessageTemplate => M_1D % Message
    
    do iM = 1, M_1D % nMessages
      M => M_1D % Message ( iM )
      call M % Initialize ( C, Tag ( iM ), Rank ( iM ), nElements ( iM ) )
    end do

    nullify ( M )

  end subroutine InitializeAllocate


  subroutine Receive ( M_1D )

    class ( MessageIncoming_1D_C_Form ), intent ( inout ), target :: &
      M_1D

    integer ( KDI ) :: &
      iM    !-- iMessage
    class ( MessageTemplate ), dimension ( : ), pointer :: &
      M

    M => M_1D % Message 
    select type ( M )
    type is ( MessageIncoming_C_Form )
      do iM = 1, M_1D % nMessages
        call M ( iM ) % Receive ( )
      end do
    end select

    nullify ( M )

  end subroutine Receive


  elemental subroutine Finalize ( M_1D )

    type ( MessageIncoming_1D_C_Form ), intent ( inout ) :: &
      M_1D

    !-- Trigger finalization of parent type

    if ( allocated ( M_1D % Message ) ) deallocate ( M_1D % Message )

  end subroutine Finalize


end module MessageIncoming_1D_C__Form
