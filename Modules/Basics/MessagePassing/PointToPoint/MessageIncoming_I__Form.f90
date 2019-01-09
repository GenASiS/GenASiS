!-- MessageIncoming_I_Form is inherited from Message_I_Form to 
!   provide specific methods for receiving messages.

module MessageIncoming_I__Form

  use MPI
  use Specifiers
  use Message_I__Form

  implicit none
  private

  type, public, extends ( Message_I_Form ) :: MessageIncoming_I_Form
  contains
    procedure, public, pass :: &
      Receive
    final :: &
      Finalize
  end type MessageIncoming_I_Form
  
contains


  subroutine Receive ( M )

    class ( MessageIncoming_I_Form ), intent ( inout ) :: &
      M

    integer :: &
      PlainInteger
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      ReceiveCount

    inquire ( iolength = ThisValueSize ) M % Value ( 1 )
    inquire ( iolength = PlainValueSize ) PlainInteger
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    ReceiveCount = size ( M % Value ) * SizeRatio
        
    call MPI_IRECV &
           ( M % Value, ReceiveCount, MPI_INTEGER, M % Rank, M % Tag, &
             M % Communicator % Handle, M % Handle, M % Error )

  end subroutine Receive


  elemental subroutine Finalize ( M )

    type ( MessageIncoming_I_Form ), intent ( inout ) :: &
      M 

    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageIncoming_I__Form

