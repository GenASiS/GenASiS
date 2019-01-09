!-- CollectiveOperation_R_Form provides a concrete extension of 
!   CollectiveOperationTemplate for real datatype for handling collective
!   operations.

module CollectiveOperation_R__Form

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
    CollectiveOperation_R_Form
      type ( MessageIncoming_R_Form ), allocatable :: &
        Incoming
      type ( MessageOutgoing_R_Form ), allocatable :: &
        Outgoing
  contains
    procedure, public, pass :: &
      InitializeAllocate
    procedure, public, pass :: &
      InitializeAssociate
    generic :: &
      Initialize => InitializeAllocate, InitializeAssociate
    procedure, public, pass :: &
      Broadcast
    procedure, public, pass :: &
      Gather
    procedure, public, pass :: &
      AllToAll
    procedure, public, pass :: &
      AllToAll_V
    procedure, public, pass :: &
      Reduce
    final :: &
      Finalize
  end type CollectiveOperation_R_Form
  
contains


  subroutine InitializeAllocate &
               ( CO, C, nOutgoing, nIncoming, RootOption )

    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
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
  
  
  subroutine InitializeAssociate &
               ( CO, C, OutgoingValue, IncomingValue, RootOption )

    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      OutgoingValue, &
      IncomingValue
    integer ( KDI ), intent ( in ), optional :: &
      RootOption

    call CO % InitializeTemplate &
           ( C, shape ( OutgoingValue ), shape ( IncomingValue ), RootOption )
    
    allocate ( CO % Incoming )
    call CO % Incoming % Initialize ( C, IncomingValue, 0, C % Rank )
    
    allocate ( CO % Outgoing )
    call CO % Outgoing % Initialize ( C, OutgoingValue, 0, C % Rank )

  end subroutine InitializeAssociate
  
  
  subroutine Broadcast ( CO )

    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO

    real :: &
      PlainReal
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
    
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( OV ) * SizeRatio
    
    call MPI_BCAST &
           ( OV, SendCount, MPI_REAL, &
             CO % Root, CO % Communicator % Handle, CO % Error)
    call Copy ( OV, IV )
            
    end associate

  end subroutine Broadcast
      

  subroutine Gather ( CO )
  
    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO
      
    real :: &
      PlainReal
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
    
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( OV ) * SizeRatio
    
    if ( CO % Root /= UNSET ) then
      call MPI_GATHER &
             ( OV, SendCount, MPI_REAL, &
               IV, SendCount, MPI_REAL, &
               CO % Root, CO % Communicator % Handle, CO % Error)  
    else
      call MPI_ALLGATHER &
             ( OV, SendCount, MPI_REAL, &
               IV, SendCount, MPI_REAL, &
               CO % Communicator % Handle, CO % Error)  
    end if

    end associate

  end subroutine Gather


  subroutine AllToAll ( CO )
  
    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO

    real :: &
      PlainReal
    integer ( KDI ) :: &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio, &
      SendCount

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
    
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = ( size ( OV ) / CO % Communicator % Size ) * SizeRatio
    
    call MPI_ALLTOALL &
           ( OV, SendCount, MPI_REAL, &
             IV, SendCount, MPI_REAL, &
             CO % Communicator % Handle, CO % Error)  
      
    end associate
    
  end subroutine AllToAll


  subroutine AllToAll_V ( CO )
  
    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO

    real :: &
      PlainReal
    integer ( KDI ) :: &
      i, &
      PlainValueSize, &
      ThisValueSize, &
      SizeRatio
    integer ( KDI ), dimension ( CO % Communicator % Size ) :: &
      SendCount, ReceiveCount, &
      SendDisplacement, ReceiveDisplacement

    associate &
      ( OV => CO % Outgoing % Value ( : ), &
        IV => CO % Incoming % Value ( : ) )
    
    inquire ( iolength = ThisValueSize ) OV ( 1 )
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )

    SendCount    = CO % nOutgoing * SizeRatio
    ReceiveCount = CO % nIncoming * SizeRatio

    SendDisplacement ( 1 )    = 0
    ReceiveDisplacement ( 1 ) = 0
    do i = 2, CO % Communicator % Size
      SendDisplacement ( i ) &
        = SendDisplacement ( i - 1 ) + SendCount ( i - 1 )
      ReceiveDisplacement ( i ) &
        = ReceiveDisplacement ( i - 1 ) + ReceiveCount ( i - 1 )
    end do

    call MPI_ALLTOALLV &
           ( OV, SendCount, SendDisplacement, MPI_REAL, &
             IV, ReceiveCount, ReceiveDisplacement, MPI_REAL, &
             CO % Communicator % Handle, CO % Error )  
      
    end associate
    
  end subroutine AllToAll_V


  subroutine Reduce ( CO, Operation )
  
    class ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO
    integer ( KDI ), intent ( in ) :: &
      Operation

    real :: &
      PlainReal
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
    inquire ( iolength = PlainValueSize ) PlainReal
    SizeRatio = max ( 1, ThisValueSize / PlainValueSize )
    SendCount = size ( OV )
    
    select case ( SizeRatio )
    case ( 1 )
      MPI_Datatype = MPI_REAL
    case ( 2 ) 
      MPI_Datatype = MPI_DOUBLE_PRECISION
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


  impure elemental subroutine Finalize ( CO )

    type ( CollectiveOperation_R_Form ), intent ( inout ) :: &
      CO

    if ( allocated ( CO % Outgoing ) ) deallocate ( CO % Outgoing )
    if ( allocated ( CO % Incoming ) ) deallocate ( CO % Incoming )
    
    if ( allocated ( CO % nOutgoing ) ) deallocate ( CO % nOutgoing )
    if ( allocated ( CO % nIncoming ) ) deallocate ( CO % nIncoming )
    
    nullify ( CO % Communicator )

  end subroutine Finalize

  
end module CollectiveOperation_R__Form
