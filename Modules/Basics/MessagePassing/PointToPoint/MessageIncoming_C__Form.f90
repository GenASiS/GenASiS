!-- MessageIncoming_C_Form is inherited from Message_C_Form to 
!   provide specific methods for receiving messages.

module MessageIncoming_C__Form

  use MPI
  use Specifiers
  use Message_C__Form

  implicit none
  private

  type, public, extends ( Message_C_Form ) :: MessageIncoming_C_Form
  contains
    procedure, public, pass :: &
      Receive
   final :: &
     Finalize
  end type MessageIncoming_C_Form

contains


  subroutine Receive ( M )

    class ( MessageIncoming_C_Form ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      ReceiveCount
    complex :: &
      PlainComplex

    inquire ( iolength = ThisValueSize ) M % Value ( 1 )
    inquire ( iolength = PlainValueSize ) PlainComplex
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    ReceiveCount = size ( M % Value ) * SizeRatio
        
    call MPI_IRECV &
           ( M % Value, ReceiveCount, MPI_COMPLEX, M % Rank, M % Tag, &
             M % Communicator % Handle, M % Handle, M % Error )

  end subroutine Receive


  elemental subroutine Finalize ( M )

    type ( MessageIncoming_C_Form ), intent ( inout ) :: &
      M 
      
    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageIncoming_C__Form
