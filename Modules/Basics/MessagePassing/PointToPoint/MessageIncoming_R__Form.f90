!-- MessageIncoming_R_Form is inherited from Message_R_Form to 
!   provide specific methods for receiving messages.

module MessageIncoming_R__Form

  use iso_c_binding
  use MPI
  use Specifiers
  use Message_R__Form

  implicit none
  private

  type, public, extends ( Message_R_Form ) :: MessageIncoming_R_Form
  contains
    procedure, public, pass :: &
      Receive
    final :: &
      Finalize
  end type MessageIncoming_R_Form

contains


  subroutine Receive ( M )

    class ( MessageIncoming_R_Form ), intent ( inout ), target :: &
      M

    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      ReceiveCount
    real :: &
      PlainReal
    real ( KDR ), dimension ( : ), pointer :: &
      Value
  
    inquire ( iolength = ThisValueSize ) M % Value ( 1 )
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    ReceiveCount = size ( M % Value ) * SizeRatio
    
    if ( c_associated ( M % D_Value ) ) then
      call c_f_pointer ( M % D_Value, Value, [ size ( M % Value ) ] )
    else
      Value => M % Value
    end if
    
    call MPI_IRECV &
           ( Value, ReceiveCount, MPI_REAL, M % Rank, M % Tag, &
             M % Communicator % Handle, M % Handle, M % Error )
    
  end subroutine Receive


  elemental subroutine Finalize ( M )

    type ( MessageIncoming_R_Form ), intent ( inout ) :: &
      M 
      
    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageIncoming_R__Form
