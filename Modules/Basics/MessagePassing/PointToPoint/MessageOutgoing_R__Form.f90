!-- MessageOutgoing_R_Form is inherited from Message_R_Form to 
!   provide specific methods for sending messages.

module MessageOutgoing_R__Form

  use MPI
  use VariableManagement
  use Message_R__Form

  implicit none
  private

  type, public, extends ( Message_R_Form ) :: MessageOutgoing_R_Form
  contains
    procedure, public, pass :: &
      Send
    final :: &
      Finalize
  end type MessageOutgoing_R_Form

contains


  subroutine Send ( M )

    class ( MessageOutgoing_R_Form ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount
    real :: &
      PlainReal
  
    inquire ( iolength = ThisValueSize ) M % Value ( 1 )
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( M % Value ) * SizeRatio

    call MPI_ISEND &
           ( M % Value, SendCount, MPI_REAL, M % Rank, M % Tag, &
             M % Communicator % Handle, M % Handle, M % Error )
    
  end subroutine Send


  elemental subroutine Finalize ( M )

    type ( MessageOutgoing_R_Form ), intent ( inout ) :: &
      M 
      
    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageOutgoing_R__Form
