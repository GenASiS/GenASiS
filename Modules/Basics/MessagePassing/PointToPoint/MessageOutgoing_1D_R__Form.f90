!-- MessageOutgoing_1D_R_Form provides the concrete extension of 
!   Message_1D_Template for real datatype to handle sending array of 
!   messages.

module MessageOutgoing_1D_R__Form

  use MPI
  use VariableManagement
  use MessagePassingBasics
  use Message_Template
  use MessageOutgoing_R__Form
  use Message_1D__Template 
 
  implicit none
  private

  type, public, extends ( Message_1D_Template ) :: MessageOutgoing_1D_R_Form
    type ( MessageOutgoing_R_Form ), dimension ( : ), allocatable :: &
      Message
  contains
    procedure, public, pass :: &
      Initialize => InitializeAllocate
    procedure, public, pass :: &
      SendOne
    procedure, public, pass :: &
      SendAll
    generic :: &
      Send => SendOne, SendAll
    final :: &
      Finalize
  end type MessageOutgoing_1D_R_Form

contains


  subroutine InitializeAllocate ( M_1D, C, Tag, Rank, nElements )

    class ( MessageOutgoing_1D_R_Form ), intent ( inout ), target :: &
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
    type ( MessageOutgoing_R_Form ), pointer :: &
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


  subroutine SendOne ( M_1D, iM )

    class ( MessageOutgoing_1D_R_Form ), intent ( inout ), target :: &
      M_1D
    integer ( KDI ), intent ( in ) :: &
      iM    !-- iMessage
    class ( MessageTemplate ), pointer :: &
      M
        
    M => M_1D % Message ( iM )
    select type ( M )
    type is ( MessageOutgoing_R_Form )
      call M % Send ( )
    end select

    nullify ( M )

  end subroutine SendOne

 
  subroutine SendAll ( M_1D )

    class ( MessageOutgoing_1D_R_Form ), intent ( inout ), target :: &
      M_1D

    integer ( KDI ) :: &
      iM    !-- iMessage
    class ( MessageTemplate ), dimension ( : ), pointer :: &
      M

    M => M_1D % Message
    select type ( M )
    type is ( MessageOutgoing_R_Form )
      do iM = 1, M_1D % nMessages
        call M ( iM ) % Send ( )
      end do
    end select

    nullify ( M )

  end subroutine SendAll
  
  
  elemental subroutine Finalize ( M_1D )
  
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      M_1D

    !-- Triggers finalization of parent type
  
    if ( allocated ( M_1D % Message ) ) deallocate ( M_1D % Message )

  end subroutine Finalize


end module MessageOutgoing_1D_R__Form
