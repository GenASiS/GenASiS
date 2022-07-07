program MessageIncomingOutgoing_Forms_Test

  use MPI
  use Specifiers
  use Display
  use Devices
  use MessagePassingBasics
  use MessageIncoming_I__Form
  use MessageIncoming_BI__Form
  use MessageIncoming_R__Form
  use MessageIncoming_C__Form
  use MessageOutgoing_I__Form
  use MessageOutgoing_BI__Form
  use MessageOutgoing_R__Form
  use MessageOutgoing_C__Form

  implicit none

  integer ( KDI ), parameter :: &
    nReceive = 8, &
    nSend = 8, &
    nReceiveBig = 1000000, &
    nSendBig = 1000000, &
    nIterationsBig = 10, &
    Tag = 100
  integer ( KDI ) :: &
    iIB, &  !-- iIterationsBig
    iV, &  !-- iValue
    SourceRank, &
    TargetRank
  real ( KDR ), dimension ( : ), allocatable :: &
    SendBuffer, &
    ReceiveBuffer
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( MessageIncoming_I_Form ), allocatable :: &
    IM_I
  type ( MessageOutgoing_I_Form ), allocatable :: &
    OM_I
  type ( MessageIncoming_BI_Form ), allocatable :: &
    IM_BI
  type ( MessageOutgoing_BI_Form ), allocatable :: &
    OM_BI
  type ( MessageIncoming_R_Form ), allocatable :: &
    IM_R
  type ( MessageOutgoing_R_Form ), allocatable :: &
    OM_R
  type ( MessageIncoming_C_Form ), allocatable :: &
    IM_C
  type ( MessageOutgoing_C_Form ), allocatable :: &
    OM_C

  associate ( PI => CONSTANT % PI )

  allocate ( C )
  
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetDisplayRank ( 0 )
  
  !-- Send to next rank, receive from previous rank
  TargetRank = modulo ( C % Rank + 1, C % Size )
  SourceRank = modulo ( C % Rank - 1, C % Size )
  
  !-- Send & Receive Integer
  
  allocate ( IM_I )
  call IM_I % Initialize ( C, Tag, SourceRank, nReceive )
  call IM_I % Receive ( ) 

  allocate ( OM_I )
  call OM_I % Initialize ( C, Tag, TargetRank, nSend )

  if ( C % Rank == 0 ) then
    OM_I % Value = [ 0, 1, 2, 3, 4, 5, 6, 7 ]
    call OM_I % Send ( )
  else
    call IM_I % Wait ( )
    OM_I % Value = IM_I % Value + 1
    call OM_I % Send ( )
  end if
 
  call Show ( OM_I % Value, 'Integer Value sent' )
  call OM_I % Wait ( )
  
  if ( C % Rank == 0 ) call IM_I % Wait ( )
  call Show ( IM_I % Value, 'Integer Value received' ) 

  deallocate ( OM_I )
  deallocate ( IM_I )

  
  !-- Send & Receive Big Integer
  
  allocate ( IM_BI )
  call IM_BI % Initialize ( C, Tag, SourceRank, nReceive )
  call IM_BI % Receive ( ) 

  allocate ( OM_BI )
  call OM_BI % Initialize ( C, Tag, TargetRank, nSend )

  if ( C % Rank == 0 ) then
    OM_BI % Value = 1_KBI * [ 0, 1, 2, 3, 4, 5, 6, 7 ]
    call OM_BI % Send ( )
  else
    call IM_BI % Wait ( )
    OM_BI % Value = IM_BI % Value + 1
    call OM_BI % Send ( )
  end if
 
  call Show ( OM_BI % Value, 'Integer ( KBI ) Value sent' )
  call OM_BI % Wait ( )
  
  if ( C % Rank == 0 ) call IM_BI % Wait ( )
  call Show ( IM_BI % Value, 'Integer ( KBI ) Value received' ) 
  
  deallocate ( OM_BI )
  deallocate ( IM_BI )

  
  !-- Send & Receive Real
  
  allocate ( IM_R )
  call IM_R % Initialize ( C, Tag, SourceRank, nReceive )
  call IM_R % Receive ( ) 

  allocate ( OM_R )
  call OM_R % Initialize ( C, Tag, TargetRank, nSend )

  if ( C % Rank == 0 ) then
    OM_R % Value = PI * [ 0, 1, 2, 3, 4, 5, 6, 7 ]
    call OM_R % Send ( )
  else
    call IM_R % Wait ( )
    OM_R % Value = IM_R % Value + PI
    call OM_R % Send ( )
  end if
 
  call Show ( OM_R % Value, 'Real Value sent' )
  call OM_R % Wait ( )
  
  if ( C % Rank == 0 ) call IM_R % Wait ( )
  call Show ( IM_R % Value, 'Real Value received' ) 
  
  deallocate ( OM_R )
  deallocate ( IM_R )

  
  !-- Send & Receive Real with associated value
  
  allocate ( ReceiveBuffer ( nReceive ) )
  allocate ( SendBuffer ( nSend ) )
  
  allocate ( IM_R )
  call IM_R % Initialize ( C, ReceiveBuffer, Tag, SourceRank )
  call IM_R % Receive ( ) 

  allocate ( OM_R )
  call OM_R % Initialize ( C, SendBuffer, Tag, TargetRank )

  if ( C % Rank == 0 ) then
    SendBuffer = PI * [ 0, 1, 2, 3, 4, 5, 6, 7 ]
    call OM_R % Send ( )
  else
    call IM_R % Wait ( )
    SendBuffer = ReceiveBuffer + PI
    call OM_R % Send ( )
  end if
 
  call Show ( SendBuffer, 'SendBuffer' )
  call OM_R % Wait ( )
  
  if ( C % Rank == 0 ) call IM_R % Wait ( )
  call Show ( ReceiveBuffer, 'ReceiveBuffer' ) 
  
  deallocate ( OM_R )
  deallocate ( IM_R )
  
  deallocate ( SendBuffer ) 
  deallocate ( ReceiveBuffer )

  !-- Send & Receive Complex
  
  allocate ( IM_C )
  call IM_C % Initialize ( C, Tag, SourceRank, nReceive )
  call IM_C % Receive ( ) 

  allocate ( OM_C )
  call OM_C % Initialize ( C, Tag, TargetRank, nSend )

  if ( C % Rank == 0 ) then
    OM_C % Value &
      = cmplx ( PI, exp ( 1.0_KDR ), kind = KDC ) * [ 0, 1, 2, 3, 4, 5, 6, 7 ]
    call OM_C % Send ( )
  else
    call IM_C % Wait ( )
    OM_C % Value = IM_C % Value + cmplx ( PI, exp ( 1.0_KDR ), kind = KDC )
    call OM_C % Send ( )
  end if
 
  call Show ( OM_C % Value, 'Complex Value sent' )
  call OM_C % Wait ( )
  
  if ( C % Rank == 0 ) call IM_C % Wait ( )
  call Show ( IM_C % Value, 'Complex Value received' ) 
  
  deallocate ( OM_C )
  deallocate ( IM_C )

  !-- Send & Receive Real, lots of big messages
  
  do iIB = 1, nIterationsBig

    allocate ( IM_R )
    call IM_R % Initialize ( C, Tag, SourceRank, nReceiveBig )
    call IM_R % AllocateDevice ( )
    call IM_R % Receive ( ) 

    allocate ( OM_R )
    call OM_R % Initialize ( C, Tag, TargetRank, nSendBig )
    call OM_R % AllocateDevice ( )

    if ( mod ( iIB, 1 ) == 0 ) then
      call Show ( iIB, 'iIB' )
      call Show ( size ( OM_R % Value ), 'nSendBig' )
      call Show ( size ( IM_R % Value ), 'nReceiveBig' )
    end if

    if ( C % Rank == 0 ) then
      OM_R % Value = PI * [ ( iV - 1, iV = 1, nSendBig ) ]
      call UpdateDevice ( OM_R % Value, OM_R % D_Value )
      call OM_R % Send ( )
    else
      call IM_R % Wait ( )
      OM_R % Value = IM_R % Value + PI
      call UpdateDevice ( OM_R % Value, OM_R % D_Value )
      call OM_R % Send ( )
    end if
 
!  call Show ( OM_R % Value, 'Real Value sent' )
    call OM_R % Wait ( )
  
    if ( C % Rank == 0 ) call IM_R % Wait ( )
!  call Show ( IM_R % Value, 'Real Value received' ) 
  
    deallocate ( OM_R )
    deallocate ( IM_R )

  end do !-- iIB

  !-- Cleanup

  deallocate ( C )

  end associate !-- PI

end program MessageIncomingOutgoing_Forms_Test
