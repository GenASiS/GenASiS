program CollectiveOperation_Form_Test

  use Specifiers
  use Display
  use MessagePassingBasics
  use REDUCTION_Singleton
  use CollectiveOperation_I__Form
  use CollectiveOperation_BI__Form
  use CollectiveOperation_R__Form
  use CollectiveOperation_C__Form

  implicit none

  integer ( KDI ) :: &
    iV, jV, &  !-- iValue, etc.
    iB, &      !-- iBuffer
    nOutgoing
  integer ( KDI ), dimension ( : ), allocatable :: &
    nOutgoing_V, &
    nIncoming_V
  complex ( KDC ) :: &
    ComplexValue
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( CollectiveOperation_I_Form ), allocatable :: &
    CO_I
  type ( CollectiveOperation_BI_Form ), allocatable :: &
    CO_BI
  type ( CollectiveOperation_R_Form ), allocatable :: &
    CO_R
  type ( CollectiveOperation_C_Form ), allocatable :: &
    CO_C

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetDisplayRank ( 0 )
  
!  nOutgoing = 120
  nOutgoing = 8

  !-- Broadcast Integer

  call Show ( 'Broadcast Integer' )

  allocate ( CO_I ) 
  call CO_I % Initialize &
         ( C, nOutgoing = [ nOutgoing ], nIncoming = [ nOutgoing ], &
           RootOption = C % Size - 1 ) 
           
  if ( C % Rank == CO_I % Root ) &
    CO_I % Outgoing % Value &
      = ( C % Rank + 1 ) * [ ( iV, iV = 1, nOutgoing ) ]
             
  call CO_I % Broadcast ( )
  
  call Show ( CO_I % Incoming % Value, 'IncomingValue_1D Integer' )

  deallocate ( CO_I ) 
    
  !-- AllGather Integer
  
  call Show ( 'AllGather Integer' )

  allocate ( CO_I ) 
  call CO_I % Initialize &
         ( C, nOutgoing = [ nOutgoing ], nIncoming = [ C % Size * nOutgoing ] ) 
           
  CO_I % Outgoing % Value &
    = ( C % Rank + 1 ) * [ ( iV, iV = 1, nOutgoing ) ]

  call CO_I % Gather ( )
  
  call Show ( CO_I % Incoming % Value, 'IncomingValue_1D Integer' )

  deallocate ( CO_I ) 
  
  !-- AllGather Integer_KBI
  
  call Show ( 'AllGather Integer_KBI' )

  allocate ( CO_BI ) 
  call CO_BI % Initialize &
         ( C, nOutgoing = [ nOutgoing ], nIncoming = [ C % Size * nOutgoing ] ) 
           
  CO_BI % Outgoing % Value &
    = int ( C % Rank + 1, KBI ) * [ ( iV, iV = 1, nOutgoing ) ]
             
  call CO_BI % Gather ( )
  
  call Show ( CO_BI % Incoming % Value, 'IncomingValue_1D Integer_KBI' )

  deallocate ( CO_BI ) 
  
  !-- AllGather Real
  
  call Show ( 'AllGather Real' )

  allocate ( CO_R ) 
  call CO_R % Initialize &
         ( C, nOutgoing = [ nOutgoing ], nIncoming = [ C % Size * nOutgoing ] ) 
           
  CO_R % Outgoing % Value &
    = ( acos ( - 1.0_KDR ) * ( C % Rank + 1 ) * [ ( iV, iV = 1, nOutgoing ) ] )
             
  call CO_R % Gather ( )

  call Show ( CO_R % Incoming % Value, 'IncomingValue_1D Real' )

  deallocate ( CO_R ) 
  
  !-- Gather Real
  
  call Show ( 'Gather Real' )

  allocate ( CO_R ) 
  call CO_R % Initialize &
         ( C, nOutgoing = [ nOutgoing ], nIncoming = [ C % Size * nOutgoing ], &
           RootOption = CONSOLE % DisplayRank ) 
           
  CO_R % Outgoing % Value &
    = ( acos ( - 1.0_KDR ) * ( C % Rank + 1 ) * [ ( iV, iV = 1, nOutgoing ) ] )
             
  call CO_R % Gather ( )

  call Show ( CO_R % Incoming % Value, 'IncomingValue_1D Real' )

  deallocate ( CO_R ) 
  
  !-- AllToAll Complex
  
  call Show ( 'AllToAll Complex' )

  allocate ( CO_C ) 
  call CO_C % Initialize & 
         ( C, nOutgoing = [ nOutgoing * C % Size ], &
           nIncoming = [ nOutgoing * C % Size ] )

  ComplexValue &
    = cmplx ( ( C % Rank + 1 ) * 1.0_KDR, &
              ( C % Rank + 1 ) * acos ( - 1.0_KDR ) )
  CO_C % Outgoing % Value &
    = spread ( ComplexValue, dim = 1, ncopies = nOutgoing * C % Size )
  
  call CO_C % AllToAll ( )
  
  call Show ( CO_C % Incoming % Value, 'IncomingValue_1D Complex' )

  deallocate ( CO_C ) 

  !-- AllToAll_V Real

  call Show ( 'AllToAll_V Real' )

  allocate ( CO_R )
  allocate ( nOutgoing_V ( C % Size ), nIncoming_V ( C % Size ) )
  nOutgoing_V = [ ( iV, iV = 1, C % Size ) ]
  nIncoming_V = [ ( nOutgoing_V ( C % Rank + 1 ), iV = 1, C % Size ) ]
  call Show ( nOutgoing_V, 'nOutgoing_V' )
  call Show ( nIncoming_V, 'nIncoming_V' )
  call CO_R % Initialize ( C, nOutgoing = nOutgoing_V, nIncoming = nIncoming_V )

  iB = 0
  do iV = 1, C % Size
    do jV = 1, iV
      iB = iB + 1
      CO_R % Outgoing % Value ( iB ) = ( C % Rank + 1 )  *  10 ** ( iV - 1 )
    end do
  end do
  call Show ( CO_R % Outgoing % Value, 'OutgoingValue_1D Real' )

  call CO_R % AllToAll_V ( )

  call Show ( CO_R % Incoming % Value, 'IncomingValue_1D Real' )

  deallocate ( CO_R )

  !-- AllReduce Real

  call Show ( 'AllReduce Real' )

  allocate ( CO_R ) 
  call CO_R % Initialize & 
         ( C, nOutgoing = [ nOutgoing ], nIncoming = [ nOutgoing ] )

  CO_R % Outgoing % Value &
    = ( C % Rank + 1 ) * [ ( real ( iV, KDR ), iV = 1, nOutgoing ) ]

  call CO_R % Reduce ( REDUCTION % SUM )
 
  call Show ( CO_R % Incoming % Value, 'IncomingValue_1D Real' )

  deallocate ( CO_R )
  
  deallocate ( C )

end program CollectiveOperation_Form_Test
