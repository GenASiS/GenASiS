!-- MessageOutgoing_I_Form is inherited from Message_I_Form to 
!   provide specific methods for sending messages.

module MessageOutgoing_I__Form

  use MPI
  use VariableManagement
  use Message_I__Form

  implicit none
  private

  type, public, extends ( Message_I_Form ) :: MessageOutgoing_I_Form
  contains
    procedure, public, pass :: &
      Send
    final :: &
      Finalize
  end type MessageOutgoing_I_Form

contains


  subroutine Send ( M )

    class ( MessageOutgoing_I_Form ), intent ( inout ) :: &
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

    type ( MessageOutgoing_I_Form ), intent ( inout ) :: &
      M 

    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageOutgoing_I__Form
