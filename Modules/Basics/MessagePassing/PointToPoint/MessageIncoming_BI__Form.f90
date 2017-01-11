!-- MessageIncoming_BI_Form is inherited from Message_BI_Form to 
!   provide specific methods for receiving messages.

module MessageIncoming_BI__Form

  use MPI
  use VariableManagement
  use Message_BI__Form

  implicit none
  private

  type, public, extends ( Message_BI_Form ) :: MessageIncoming_BI_Form
  contains
    procedure, public, pass :: &
      Receive
    final :: &
      Finalize
  end type MessageIncoming_BI_Form

contains


  subroutine Receive ( M )

    class ( MessageIncoming_BI_Form ), intent ( inout ) :: &
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

    type ( MessageIncoming_BI_Form ), intent ( inout ) :: &
      M 
      
    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageIncoming_BI__Form
