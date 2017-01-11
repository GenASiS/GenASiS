program IncomingOutgoingMessageArray_Forms_Test

  use VariableManagement
  use Display
  use MessagePassing
  use FileSystem
  use PROGRAM_HEADER_Singleton

  implicit none
  
  integer ( KDI ) :: &
    iR, &  !-- iRank
    iM, &  !-- iMessage
    iV, &  !-- iValue
    iIB, &  !-- iIterationsBig
    nIterationsBig
  integer ( KDI ), dimension ( : ), allocatable :: &
    nReceive, &
    nSend, &
    SourceRank, &
    TargetRank, &
    Tag
  character ( LDF ) :: &
    Name = 'IncomingOutgoingMessageArray_Forms_Test'
  logical ( KDL ) :: &
    AllFinished
  type ( IncomingMessageArrayIntegerForm ), allocatable :: &
    IMA_I
  type ( OutgoingMessageArrayIntegerForm ), allocatable :: &
    OMA_I
  type ( IncomingMessageArrayBigIntegerForm ), allocatable :: &
    IMA_BI
  type ( OutgoingMessageArrayBigIntegerForm ), allocatable :: &
    OMA_BI
  type ( IncomingMessageArrayRealForm ), allocatable :: &
    IMA_R
  type ( OutgoingMessageArrayRealForm ), allocatable :: &
    OMA_R
  type ( IncomingMessageArrayComplexForm ), allocatable :: &
    IMA_C
  type ( OutgoingMessageArrayComplexForm ), allocatable :: &
    OMA_C
  
  !-- Recreate AllToAll functionality using Send & Receive

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize &
         ( Name, AppendDimensionalityOption = .false. )
  
  allocate ( nReceive ( PROGRAM_HEADER % Communicator % Size ) )
  allocate ( nSend ( PROGRAM_HEADER % Communicator % Size ) )
  allocate ( SourceRank ( PROGRAM_HEADER % Communicator % Size ) )
  allocate ( TargetRank ( PROGRAM_HEADER % Communicator % Size ) )
  allocate ( Tag ( PROGRAM_HEADER % Communicator % Size ) )
  
  nReceive = 4
  nSend = 4
  SourceRank = [ ( iR, iR = 0, PROGRAM_HEADER % Communicator % Size - 1 ) ]
  TargetRank = [ ( iR, iR = 0, PROGRAM_HEADER % Communicator % Size - 1 ) ]
  Tag = 100
  
  !-- test integer
  
  allocate ( IMA_I )
  allocate ( OMA_I )

  call IMA_I % Initialize ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
  call OMA_I % Initialize ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )
  
  call IMA_I % Receive ( ) 
  
  do iM = 1, OMA_I % nMessages
    OMA_I % Message ( iM ) % Value &
      = ( PROGRAM_HEADER % Communicator % Rank + 1 ) &
        * [ ( iV, iV = 1, nSend ( iM ) ) ]
  end do

  call OMA_I % Send ( )
  
  AllFinished = .false. 
  do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_I % Wait ( AllFinished, iR )
    call IMA_I % WaitAny ( AllFinished, iR )
    if ( AllFinished ) exit
    call Show ( iR, 'iR', nLeadingLinesOption = 1 )
    call Show ( IMA_I % Message ( iR ) % Value, 'ReceivedInteger' )
  end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_I % Wait ( )
  call OMA_I % WaitAll ( )
  
  deallocate ( OMA_I )
  deallocate ( IMA_I )
  
  !-- test integer KBI
  
  allocate ( IMA_BI )
  allocate ( OMA_BI )

  call IMA_BI % Initialize ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
  call OMA_BI % Initialize ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

  call IMA_BI % Receive ( ) 
  
  do iM = 1, OMA_BI % nMessages
    OMA_BI % Message ( iM ) % Value &
      = ( PROGRAM_HEADER % Communicator % Rank + 1 ) * [ ( int ( iV, KBI ), iV = 1, nSend ( iM ) ) ]
  end do

  call OMA_BI % Send ( )
  
  AllFinished = .false. 
  do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_BI % Wait ( AllFinished, iR )
    call IMA_BI % WaitAny ( AllFinished, iR )
    if ( AllFinished ) exit
    call Show ( iR, 'iR', nLeadingLinesOption = 1 )
    call Show ( IMA_BI % Message ( iR ) % Value, 'ReceivedBigInteger' )
  end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_BI % Wait ( )
  call OMA_BI % WaitAll ( )
  
  deallocate ( OMA_BI )
  deallocate ( IMA_BI )
    
  !-- test real
  
  allocate ( IMA_R )
  allocate ( OMA_R )

  call IMA_R % Initialize ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
  call OMA_R % Initialize ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

  call IMA_R % Receive ( ) 
  
  do iM = 1, OMA_R % nMessages
    OMA_R % Message ( iM ) % Value &
      = ( PROGRAM_HEADER % Communicator % Rank + 1 ) * [ ( real ( iV, KDR ), iV = 1, nSend ( iM ) ) ]
  end do

  call OMA_R % Send ( )
  
  AllFinished = .false. 
  do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_R % Wait ( AllFinished, iR )
    call IMA_R % WaitAny ( AllFinished, iR )
    if ( AllFinished ) exit
    call Show ( iR, 'iR', nLeadingLinesOption = 1 )
    call Show ( IMA_R % Message ( iR ) % Value, 'ReceivedReal' )
  end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_R % Wait ( )
  call OMA_R % WaitAll ( )
  
  deallocate ( OMA_R )
  deallocate ( IMA_R )    

  !-- test complex
  
  allocate ( IMA_C )
  allocate ( OMA_C )

  call IMA_C % Initialize ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
  call OMA_C % Initialize ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

  call IMA_C % Receive ( ) 
  
  do iM = 1, OMA_C % nMessages
    OMA_C % Message ( iM ) % Value &
      = ( PROGRAM_HEADER % Communicator % Rank + 1 ) &
        * [ ( cmplx ( iV, 2 * iV, KDC ), iV = 1, nSend ( iM ) ) ]
  end do

  call OMA_C % Send ( )
  
  AllFinished = .false. 
  do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_C % Wait ( AllFinished, iR )
    call IMA_C % WaitAny ( AllFinished, iR )
    if ( AllFinished ) exit
    call Show ( iR, 'iR', nLeadingLinesOption = 1 )
    call Show ( IMA_C % Message ( iR ) % Value, 'ReceivedComplex' )
  end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_C % Wait ( )
  call OMA_C % WaitAll ( )
  
  deallocate ( OMA_C )
  deallocate ( IMA_C )    


  !-- test lots of big messages

  nReceive = 1000000
  nSend = 1000000
  nIterationsBig = 1000
  
  do iIB = 1, nIterationsBig

    if ( mod ( iIB, 1 ) == 0 ) then
      call Show ( iIB, 'iIB' )
      call Show ( nSend, 'nSendBig' )
      call Show ( nReceive, 'nReceiveBig' )
      call PROGRAM_HEADER % ShowStatistics &
             ( CONSOLE % INFO_1, AcrossProcessesOption = .true. )
    end if

    !-- Integer

    allocate ( IMA_I )
    allocate ( OMA_I )

    call IMA_I % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
    call OMA_I % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

    call IMA_I % Receive ( ) 
  
    do iM = 1, OMA_I % nMessages
      OMA_I % Message ( iM ) % Value &
        = ( PROGRAM_HEADER % Communicator % Rank + 1 ) &
          * [ ( iV, iV = 1, nSend ( iM ) ) ]
    end do

    call OMA_I % Send ( )
  
    AllFinished = .false. 
    do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_I % Wait ( AllFinished, iR )
      call IMA_I % WaitAny ( AllFinished, iR )
      if ( AllFinished ) exit
!      call Show ( iR, 'iR', nLeadingLinesOption = 1 )
!      call Show ( IMA_I % Message ( iR ) % Value, 'ReceivedReal' )
    end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_I % Wait ( )
    call OMA_I % WaitAll ( )
  
    deallocate ( OMA_I )
    deallocate ( IMA_I )    

    !-- Big Integer

    allocate ( IMA_BI )
    allocate ( OMA_BI )

    call IMA_BI % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
    call OMA_BI % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

    call IMA_BI % Receive ( ) 
  
    do iM = 1, OMA_BI % nMessages
      OMA_BI % Message ( iM ) % Value &
        = ( PROGRAM_HEADER % Communicator % Rank + 1 ) &
          * [ ( int ( iV, KBI ), iV = 1, nSend ( iM ) ) ]
    end do

    call OMA_BI % Send ( )
  
    AllFinished = .false. 
    do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_BI % Wait ( AllFinished, iR )
      call IMA_BI % WaitAny ( AllFinished, iR )
      if ( AllFinished ) exit
!      call Show ( iR, 'iR', nLeadingLinesOption = 1 )
!      call Show ( IMA_BI % Message ( iR ) % Value, 'ReceivedReal' )
    end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_BI % Wait ( )
    call OMA_BI % WaitAll ( )
  
    deallocate ( OMA_BI )
    deallocate ( IMA_BI )    

    !-- Real

    allocate ( IMA_R )
    allocate ( OMA_R )

    call IMA_R % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
    call OMA_R % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

    call IMA_R % Receive ( ) 
  
    do iM = 1, OMA_R % nMessages
      OMA_R % Message ( iM ) % Value &
        = ( PROGRAM_HEADER % Communicator % Rank + 1 ) &
          * [ ( real ( iV, KDR ), iV = 1, nSend ( iM ) ) ]
    end do

    call OMA_R % Send ( )
  
    AllFinished = .false. 
    do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_R % Wait ( AllFinished, iR )
      call IMA_R % WaitAny ( AllFinished, iR )
      if ( AllFinished ) exit
!      call Show ( iR, 'iR', nLeadingLinesOption = 1 )
!      call Show ( IMA_R % Message ( iR ) % Value, 'ReceivedReal' )
    end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_R % Wait ( )
    call OMA_R % WaitAll ( )
  
    deallocate ( OMA_R )
    deallocate ( IMA_R )    

    !-- Complex

    allocate ( IMA_C )
    allocate ( OMA_C )

    call IMA_C % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, SourceRank, nReceive )
    call OMA_C % Initialize &
           ( PROGRAM_HEADER % Communicator, Tag, TargetRank, nSend )

    call IMA_C % Receive ( ) 
  
    do iM = 1, OMA_C % nMessages
      OMA_C % Message ( iM ) % Value &
        = ( PROGRAM_HEADER % Communicator % Rank + 1 ) &
          * [ ( cmplx ( iV, 2 * iV, KDC ), iV = 1, nSend ( iM ) ) ]
    end do

    call OMA_C % Send ( )
  
    AllFinished = .false. 
    do while ( .not. AllFinished )
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!    call IMA_C % Wait ( AllFinished, iR )
      call IMA_C % WaitAny ( AllFinished, iR )
      if ( AllFinished ) exit
!      call Show ( iR, 'iR', nLeadingLinesOption = 1 )
!      call Show ( IMA_C % Message ( iR ) % Value, 'ReceivedReal' )
    end do
  
!-- FIXME: Generic overloading not inherited by extension in Intel 12.1.3
!  call OMA_C % Wait ( )
    call OMA_C % WaitAll ( )
  
    deallocate ( OMA_C )
    deallocate ( IMA_C )    

  end do !-- iIB

  !-- Cleanup

  deallocate ( PROGRAM_HEADER )  
  
end program IncomingOutgoingMessageArray_Forms_Test
