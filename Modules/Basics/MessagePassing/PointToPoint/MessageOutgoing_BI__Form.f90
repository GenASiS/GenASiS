!-- MessageOutgoing_BI_Form is inherited from Message_BI_Form to 
!   provide specific methods for sending messages.

module MessageOutgoing_BI__Form

  use MPI
  use VariableManagement
  use Message_BI__Form

  implicit none
  private

  type, public, extends ( Message_BI_Form ) :: MessageOutgoing_BI_Form
  contains
    procedure, public, pass :: &
      Send => Send
    final :: &
      Finalize
  end type MessageOutgoing_BI_Form

contains


  subroutine Send ( M )

    class ( MessageOutgoing_BI_Form ), intent ( inout ) :: &
      M

    integer :: &
      PlainInteger
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    inquire ( iolength = ThisValueSize ) M % Value ( 1 )
    inquire ( iolength = PlainValueSize ) PlainInteger
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( M % Value ) * SizeRatio

    call MPI_ISEND &
           ( M % Value, SendCount, MPI_INTEGER, M % Rank, M % Tag, &
             M % Communicator % Handle, M % Handle, M % Error )
    
  end subroutine Send


  elemental subroutine Finalize ( M )

    type ( MessageOutgoing_BI_Form ), intent ( inout ) :: &
      M 
      
    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageOutgoing_BI__Form
