!-- MessageOutgoing_C_Form is inherited from Message_C_Form to 
!   provide specific methods for sending messages.

module MessageOutgoing_C__Form

  use MPI
  use Specifiers
  use Message_C__Form

  implicit none
  private

  type, public, extends ( Message_C_Form ) :: MessageOutgoing_C_Form
  contains
    procedure, public, pass :: &
      Send
    final :: &
      Finalize
  end type MessageOutgoing_C_Form

contains


  subroutine Send ( M )

    class ( MessageOutgoing_C_Form ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount
    complex :: &
      PlainComplex

    inquire ( iolength = ThisValueSize ) M % Value ( 1 )
    inquire ( iolength = PlainValueSize ) PlainComplex
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( M % Value ) * SizeRatio
        
    call MPI_ISEND &
           ( M % Value, SendCount, MPI_COMPLEX, M % Rank, M % Tag, &
             M % Communicator % Handle, M % Handle, M % Error )

  end subroutine Send


  elemental subroutine Finalize ( M )

    type ( MessageOutgoing_C_Form ), intent ( inout ) :: &
      M 
      
    !-- Trigger finalization of parent type
    
  end subroutine Finalize
  
  
end module MessageOutgoing_C__Form
