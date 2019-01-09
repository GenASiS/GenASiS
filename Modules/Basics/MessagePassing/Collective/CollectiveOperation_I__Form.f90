!-- CollectiveOperation_I_Form provides a concrete extension of 
!   CollectiveOperationTemplate for integer datatype for handling collective
!   operations.

module CollectiveOperation_I__Form

  use MPI
  use Specifiers
  use DataManagement
  use Display
  use MessagePassingBasics
  use PointToPoint
  use CollectiveOperation_Template

  implicit none
  private

  type, public, extends ( CollectiveOperationTemplate ) :: &
    CollectiveOperation_I_Form
      type ( MessageIncoming_I_Form ), allocatable :: &
        Incoming
      type ( MessageOutgoing_I_Form ), allocatable :: &
        Outgoing
  contains
    procedure, public, pass :: &
      Initialize => InitializeAllocate
    procedure, public, pass :: &
      Broadcast
    procedure, public, pass :: &
      Gather
    procedure, public, pass :: &
      AllToAll
    procedure, public, pass :: &
      Reduce
    final :: &
      Finalize
  end type CollectiveOperation_I_Form
  
contains


  subroutine InitializeAllocate &
               ( CO, C, nOutgoing, nIncoming, RootOption )

    class ( CollectiveOperation_I_Form ), intent ( inout ) :: &
      CO
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ), dimension ( : ) :: &
      nOutgoing, &
      nIncoming
    integer ( KDI ), intent ( in ), optional :: &
      RootOption
      
    call CO % InitializeTemplate ( C, nOutgoing, nIncoming, RootOption )

    allocate ( CO % Incoming )
    call CO % Incoming % Initialize ( C, 0, C % Rank, sum ( nIncoming ) )
    
    allocate ( CO % Outgoing )
    call CO % Outgoing % Initialize ( C, 0, C % Rank, sum ( nOutgoing ) )
      
  end subroutine InitializeAllocate
  
  
  subroutine Broadcast ( CO )

    class ( CollectiveOperation_I_Form ), intent ( inout ) :: &
      CO

    integer :: &
      PlainInteger
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
      
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainInteger
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( OV ) * SizeRatio
      
    call MPI_BCAST &
           ( OV, SendCount, MPI_INTEGER, &
             CO % Root, CO % Communicator % Handle, CO % Error)
    call Copy ( OV, IV )
            
    end associate

  end subroutine Broadcast


  subroutine Gather ( CO )
  
    class ( CollectiveOperation_I_Form ), intent ( inout ) :: &
      CO
      
    integer :: &
      PlainInteger
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
      
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainInteger
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( OV ) * SizeRatio
      
    if ( CO % Root /= UNSET ) then
      call MPI_GATHER &
             ( OV, SendCount, MPI_INTEGER, &
               IV, SendCount, MPI_INTEGER, &
               CO % Root, CO % Communicator % Handle, CO % Error)
    else
      call MPI_ALLGATHER &
             ( OV, SendCount, MPI_INTEGER, &
               IV, SendCount, MPI_INTEGER, &
               CO % Communicator % Handle, CO % Error)  
    end if
      
    end associate

  end subroutine Gather


  subroutine AllToAll ( CO )
  
    class ( CollectiveOperation_I_Form ), intent ( inout ) :: &
      CO

    integer :: &
      PlainInteger
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
      
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainInteger
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = ( size ( OV ) / CO % Communicator % Size ) * SizeRatio
      
    call MPI_ALLTOALL &
           ( OV, SendCount, MPI_INTEGER, &
             IV, SendCount, MPI_INTEGER, &
             CO % Communicator % Handle, CO % Error)  
      
    end associate
    
  end subroutine AllToAll


  subroutine Reduce ( CO, Operation )
  
    class ( CollectiveOperation_I_Form ), intent ( inout ) :: &
      CO
    integer ( KDI ), intent ( in ) :: &
      Operation

    integer :: &
      PlainInteger
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount, &
      MPI_Datatype
    
    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
      
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainInteger
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( OV )

    select case ( SizeRatio )
    case ( 1 )
      MPI_Datatype = MPI_INTEGER
    case default
      call Show &
             ( 'MPI datatype error in CollectiveOperation_Form % Reduce', &
               CONSOLE % ERROR )
    end select
      
    if ( CO % Root /= UNSET ) then
      call MPI_REDUCE &
             ( OV, IV, SendCount, MPI_Datatype, Operation, &
               CO % Root, CO % Communicator % Handle, CO % Error )
    else
      call MPI_ALLREDUCE &
             ( OV, IV, SendCount, MPI_Datatype, Operation, &
               CO % Communicator % Handle, CO % Error )
    end if
      
    end associate

  end subroutine Reduce


  elemental subroutine Finalize ( CO )

    type ( CollectiveOperation_I_Form ), intent ( inout ) :: &
      CO

    if ( allocated ( CO % Outgoing ) ) deallocate ( CO % Outgoing )
    if ( allocated ( CO % Incoming ) ) deallocate ( CO % Incoming )
    
    if ( allocated ( CO % nOutgoing ) ) deallocate ( CO % nOutgoing )
    if ( allocated ( CO % nIncoming ) ) deallocate ( CO % nIncoming )
    
    nullify ( CO % Communicator )

  end subroutine Finalize

  
end module CollectiveOperation_I__Form
