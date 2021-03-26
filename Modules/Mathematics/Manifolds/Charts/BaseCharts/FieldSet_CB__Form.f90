module FieldSet_CB__Form

  !-- FieldSet_ChartBase__Form

  use Basics
  use ManifoldBasics
  use ChartBasics
  use Chart_BH__Form

  implicit none
  private

  type, public, extends ( FieldSet_CH_Form ) :: FieldSet_CB_Form
    logical ( KDL ) :: &
      UseDeviceGhost
    integer ( KDI ) :: &
      iTimerGhostCommunication = 0, &
      iTimerGhostPackUnpack    = 0, &
      nValues  = 0, &
      nFields  = 0
    character ( LDL ), dimension ( : ), allocatable :: &
      Field
    class ( StorageForm ), allocatable :: &
      FieldSet, &
      FieldSetStream
    type ( MessageIncoming_1D_R_Form ), allocatable :: &
      IncomingFace_L_R, &
      IncomingFace_R_L, &
      IncomingEdge_LL_RR, &
      IncomingEdge_RR_LL, &
      IncomingEdge_LR_RL, &
      IncomingEdge_RL_LR
    type ( MessageOutgoing_1D_R_Form ), allocatable :: &
      OutgoingFace_L_R, &
      OutgoingFace_R_L, &
      OutgoingEdge_LL_RR, &
      OutgoingEdge_RR_LL, &
      OutgoingEdge_LR_RL, &
      OutgoingEdge_RL_LR
    type ( Stream_CH_Pointer ), dimension ( : ), allocatable :: &
      Stream
  contains
    procedure, public, pass :: &
      InitializeAllocate
    generic, public :: &
      Initialize => InitializeAllocate
    procedure, public, pass :: &
      ExchangeGhostData
    procedure, public, pass :: &
      StartGhostExchange
    procedure, public, pass :: &
      FinishGhostExchange
    procedure, public, pass :: &
      AddStream
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateFieldSet
  end type FieldSet_CB_Form

    private :: &
      StartExchangeFace, &
      FinishExchangeFace, &
      StartExchangeEdge, &
      FinishExchangeEdge

      private :: &
        LoadMessage, &
        StoreMessage

    integer ( KDI ), private, parameter :: &
      MAX_STREAMS = MANIFOLD % MAX_STREAMS
    integer ( KDI ), dimension ( 3 ), private, parameter :: &
      !-- Faces
      TAG_RECEIVE_FACE_L  = [ 99, 98, 97 ], &
      TAG_RECEIVE_FACE_R  = [ 96, 95, 94 ], &
      TAG_SEND_FACE_L     = TAG_RECEIVE_FACE_R, &
      TAG_SEND_FACE_R     = TAG_RECEIVE_FACE_L, &
      !-- Edges
      TAG_RECEIVE_EDGE_LL = [ 93, 92, 91 ], &
      TAG_RECEIVE_EDGE_RR = [ 90, 89, 88 ], &
      TAG_RECEIVE_EDGE_LR = [ 87, 86, 85 ], &
      TAG_RECEIVE_EDGE_RL = [ 84, 83, 82 ], &
      TAG_SEND_EDGE_LL    = TAG_RECEIVE_EDGE_RR, &
      TAG_SEND_EDGE_RR    = TAG_RECEIVE_EDGE_LL, &
      TAG_SEND_EDGE_LR    = TAG_RECEIVE_EDGE_RL, &
      TAG_SEND_EDGE_RL    = TAG_RECEIVE_EDGE_LR


contains


  subroutine InitializeAllocate &
               ( FSC, C, FSM, NameShort, nFields, FieldOption, PinnedOption, &
                 UseDeviceGhostOption, IgnorabilityOption )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_BH_Form ), intent ( inout ) :: &
      C
    class ( FieldSet_MH_Form ), intent ( in ) :: &
      FSM
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    logical ( KDL ), intent ( in ), optional :: &
      PinnedOption, &
      UseDeviceGhostOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FSC % Type == '' ) &
      FSC % Type = 'a FieldSet_CB' 

    FSC % nValues  =  C % nValues
    FSC % nFields  =  nFields

    allocate ( FSC % Field ( nFields ) )
    FSC % Field = ''
    if ( present ( FieldOption ) ) &
      FSC % Field = FieldOption   

    call FSC % Initialize ( C, FSM )

    call FSC % AllocateFieldSet ( )

    FSC % UseDeviceGhost  =  .false.
    if ( present ( UseDeviceGhostOption ) )  &
      FSC % UseDeviceGhost  =  UseDeviceGhostOption  

    allocate ( FSC % Stream ( MAX_STREAMS ) )

  end subroutine InitializeAllocate


  subroutine ExchangeGhostData  ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    call FSC % StartGhostExchange ( )
    call FSC % FinishGhostExchange ( )
   
  end subroutine ExchangeGhostData


  subroutine StartGhostExchange  ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    !-- Start faces
    call StartExchangeFace &
           ( FSC, FSC % IncomingFace_L_R, FSC % OutgoingFace_L_R, &
             C % PortalFace_L_R, TAG_RECEIVE_FACE_L, TAG_SEND_FACE_R )
    call StartExchangeFace &
           ( FSC, FSC % IncomingFace_R_L, FSC % OutgoingFace_R_L, &
             C % PortalFace_R_L, TAG_RECEIVE_FACE_R, TAG_SEND_FACE_L )

    !-- Start edges
    call StartExchangeEdge &
           ( FSC, FSC % IncomingEdge_LL_RR, FSC % OutgoingEdge_LL_RR, &
             C % PortalEdge_LL_RR, TAG_RECEIVE_EDGE_LL, TAG_SEND_EDGE_RR )
    call StartExchangeEdge &
           ( FSC, FSC % IncomingEdge_RR_LL, FSC % OutgoingEdge_RR_LL, &
             C % PortalEdge_RR_LL, TAG_RECEIVE_EDGE_RR, TAG_SEND_EDGE_LL )
    call StartExchangeEdge &
           ( FSC, FSC % IncomingEdge_LR_RL, FSC % OutgoingEdge_LR_RL, &
             C % PortalEdge_LR_RL, TAG_RECEIVE_EDGE_LR, TAG_SEND_EDGE_RL )
    call StartExchangeEdge &
           ( FSC, FSC % IncomingEdge_RL_LR, FSC % OutgoingEdge_RL_LR, &
             C % PortalEdge_RL_LR, TAG_RECEIVE_EDGE_RL, TAG_SEND_EDGE_LR )

    end select  !-- C

  end subroutine StartGhostExchange


  subroutine FinishGhostExchange  ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

!  subroutine FinishExchangeFace &
!               ( FSC, IncomingFace, OutgoingFace, TagReceive )
    !-- Finish faces
    call FinishExchangeFace &
           ( FSC, FSC % IncomingFace_L_R, FSC % OutgoingFace_L_R, &
             TAG_RECEIVE_FACE_L )
    call FinishExchangeFace &
           ( FSC, FSC % IncomingFace_R_L, FSC % OutgoingFace_R_L, &
             TAG_RECEIVE_FACE_R )

    !-- Finish edges
    call FinishExchangeEdge &
           ( FSC, FSC % IncomingEdge_LL_RR, FSC % OutgoingEdge_LL_RR, &
             TAG_RECEIVE_EDGE_LL )
    call FinishExchangeEdge &
           ( FSC, FSC % IncomingEdge_RR_LL, FSC % OutgoingEdge_RR_LL, &
             TAG_RECEIVE_EDGE_RR )
    call FinishExchangeEdge &
           ( FSC, FSC % IncomingEdge_LR_RL, FSC % OutgoingEdge_LR_RL, &
             TAG_RECEIVE_EDGE_LR )
    call FinishExchangeEdge &
           ( FSC, FSC % IncomingEdge_RL_LR, FSC % OutgoingEdge_RL_LR, &
             TAG_RECEIVE_EDGE_RL )
    
  end subroutine FinishGhostExchange


  subroutine AddStream ( FSC, SC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    class ( Stream_CH_Form ), intent ( in ), target :: &
      SC
    
    integer ( KDI ) :: &
      iS

    associate ( nS  =>  FSC % nStreams )

    do iS  =  1, nS
      if ( associated ( FSC % Stream ( iS ) % Pointer, SC ) ) then
        call Show ( 'Stream already added to ' // FSC % Type, &
                    CONSOLE % WARNING )
        call Show ( FSC % Name, 'FieldSet', CONSOLE % WARNING )
        call Show (  SC % Name, 'Stream',   CONSOLE % WARNING )
        return
      end if
    end do !-- iS

    nS  =  nS + 1
    FSC % Stream ( iS ) % Pointer  =>  SC
    call Show ( 'Adding a Stream to ' // trim ( FSC % Type ), &
                FSC % IGNORABILITY + 1 )
    call Show ( FSC % Name, 'FieldSet', FSC % IGNORABILITY + 1 )
    call Show (  SC % Name, 'Stream',   FSC % IGNORABILITY + 1 )

    end associate !-- nS

  end subroutine AddStream


  impure elemental subroutine Finalize ( FSC )

    type ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    if ( allocated ( FSC % Stream ) ) &
      deallocate ( FSC % Stream )

    if ( allocated ( FSC % OutgoingEdge_RL_LR ) ) &
      deallocate ( FSC % OutgoingEdge_RL_LR )
    if ( allocated ( FSC % OutgoingEdge_LR_RL ) ) &
      deallocate ( FSC % OutgoingEdge_LR_RL )
    if ( allocated ( FSC % OutgoingEdge_RR_LL ) ) &
      deallocate ( FSC % OutgoingEdge_RR_LL )
    if ( allocated ( FSC % OutgoingEdge_LL_RR ) ) &
      deallocate ( FSC % OutgoingEdge_LL_RR )
    if ( allocated ( FSC % OutgoingFace_R_L ) ) &
      deallocate ( FSC % OutgoingFace_R_L )
    if ( allocated ( FSC % OutgoingFace_L_R ) ) &
      deallocate ( FSC % OutgoingFace_L_R )

    if ( allocated ( FSC % IncomingEdge_RL_LR ) ) &
      deallocate ( FSC % IncomingEdge_RL_LR )
    if ( allocated ( FSC % IncomingEdge_LR_RL ) ) &
      deallocate ( FSC % IncomingEdge_LR_RL )
    if ( allocated ( FSC % IncomingEdge_RR_LL ) ) &
      deallocate ( FSC % IncomingEdge_RR_LL )
    if ( allocated ( FSC % IncomingEdge_LL_RR ) ) &
      deallocate ( FSC % IncomingEdge_LL_RR )
    if ( allocated ( FSC % IncomingFace_R_L ) ) &
      deallocate ( FSC % IncomingFace_R_L )
    if ( allocated ( FSC % IncomingFace_L_R ) ) &
      deallocate ( FSC % IncomingFace_L_R )

    if ( allocated ( FSC % FieldSetStream ) ) &
      deallocate ( FSC % FieldSetStream )
    if ( allocated ( FSC % FieldSet ) ) &
      deallocate ( FSC % FieldSet )

    if ( allocated ( FSC % Field ) ) &
      deallocate ( FSC % Field )

  end subroutine Finalize


  subroutine AllocateFieldSet ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    if ( allocated ( FSC % FieldSet ) ) &
      return

    call Show ( 'Allocating a FieldSet',    FSC % IGNORABILITY + 1 )
    call Show ( FSC % Name,    'Name',    FSC % IGNORABILITY + 1 )
    call Show ( FSC % Field,   'Field',   FSC % IGNORABILITY + 1 )
    call Show ( FSC % nFields, 'nFields', FSC % IGNORABILITY + 1 )
    call Show ( FSC % nValues, 'nValues', FSC % IGNORABILITY + 1 )
    
    allocate ( FSC % FieldSet )
    call FSC % FieldSet % Initialize &
           ( [ FSC % nValues, FSC % nFields ], &
             VariableOption = FSC % Field, NameOption = FSC % Name, &
             PinnedOption = FSC % Pinned )

    allocate ( FSC % FieldSetStream )
    call FSC % FieldSetStream % Initialize ( FSC % FieldSet )

  end subroutine AllocateFieldSet


  subroutine StartExchangeFace &
               ( FSC, IncomingFace, OutgoingFace, PH, TagReceive, TagSend )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingFace
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend
    
    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend
    type ( TimerForm ), pointer :: &
      T 
      
    T  =>  PROGRAM_HEADER % TimerPointer ( FSC % iTimerGhostCommunication )
    
    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate &
      ( Communicator  =>  C % Manifold % Communicator, &
        nCB  =>  C % nCellsBrick, &
        nGL  =>  C % nGhostLayers, &
        nD   =>  C % nDimensions )

    !-- Allocate on first use

    if ( .not. allocated ( IncomingFace ) &
         .and. .not. allocated ( OutgoingFace ) ) then
    
      allocate ( IncomingFace )
      allocate ( OutgoingFace )

      call IncomingFace % Initialize &
             ( Communicator, TagReceive ( : nD ), PH % Source, &
               PH % nChunksFrom  *  FSC % nFields )
      call OutgoingFace % Initialize &
             ( Communicator, TagSend ( : nD ), PH % Target, &
               PH % nChunksTo  *  FSC % nFields )
    
      if ( FSC % UseDeviceGhost ) then
        call IncomingFace % AllocateDevice ( )
        call OutgoingFace % AllocateDevice ( )
      end if 
    
    end if  !-- allocated faces
    
    !-- Post Receives
    
    call T % Start ( )
    call IncomingFace % Receive ( )
    call T % Stop ( )

    !-- Post Sends

    do iD = 1, nD

      nSend         =  nCB
      nSend ( iD )  =  nGL ( iD )

      !-- In setting oSend, note Copy command does not inherit lbound
      if ( TagSend ( iD )  ==  TAG_SEND_FACE_R ( iD ) ) then
        oSend         =  nGL
        oSend ( iD )  =  oSend ( iD )  +  nCB ( iD )  -  nGL ( iD )
      else if ( TagSend ( iD )  ==  TAG_SEND_FACE_L ( iD ) ) then
        oSend  =  nGL
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_CB__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartExchangeFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      call LoadMessage &
             ( FSC, OutgoingFace % Message ( iD ), nSend, oSend )
      
      call T % Start ( )
      call OutgoingFace % Send ( iD )
      call T % Stop ( )

    end do !-- iD

    !-- Cleanup

    end associate  !-- Communicator, etc.
    end select  !-- C

    nullify ( T )

  end subroutine StartExchangeFace


  subroutine FinishExchangeFace &
               ( FSC, IncomingFace, OutgoingFace, TagReceive )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    type ( MessageIncoming_1D_R_Form ), intent ( inout ) :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      OutgoingFace
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive

    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oReceive, &
      nReceive
    logical ( KDL ) :: &
      AllFinished
    type ( TimerForm ), pointer :: &
      T 
      
    T => PROGRAM_HEADER % TimerPointer ( FSC % iTimerGhostCommunication )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate &
      ( nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        iaB => C % iaBrick, &
         nB => C % nBricks )

    !-- Wait for Receives

    do 

      call T % Start ( )
      call IncomingFace % Wait ( AllFinished, iD )
      call T % Stop ( )
      
      if ( AllFinished ) exit

      nReceive        = nCB
      nReceive ( iD ) = nGL ( iD )

      !-- In setting oReceive, note Copy command does not inherit lbound
      if ( TagReceive ( iD )  ==  TAG_RECEIVE_FACE_L ( iD ) ) then
        if ( iaB ( iD )  ==  1 .and. .not. C % Periodic ( iD ) ) &
          cycle
        oReceive        =  nGL
        oReceive ( iD ) =  oReceive ( iD )  -  nGL ( iD )
      else if ( TagReceive ( iD )  ==  TAG_RECEIVE_FACE_R ( iD ) ) then
        if ( iaB ( iD )  ==  nB ( iD ) .and. .not. C % Periodic ( iD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_CB__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchangeFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage &
             ( FSC, IncomingFace % Message ( iD ), nReceive, oReceive )

    end do

    !-- Wait for Sends
    call T % Start ( )
    call OutgoingFace % Wait ( )
    call T % Stop ( )

    !-- Cleanup

    end associate !-- nCB etc.
    end select  !-- C
    
    nullify ( T )

  end subroutine FinishExchangeFace


  subroutine StartExchangeEdge &
               ( FSC, IncomingEdge, OutgoingEdge, PH, TagReceive, TagSend )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingEdge
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend
    
    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension, etc.
      kM  !-- kMessage
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend
    logical ( KDL ), dimension ( 3 ) :: &
      DimensionMask
    type ( TimerForm ), pointer :: &
      T 
      
    T => PROGRAM_HEADER % TimerPointer ( FSC % iTimerGhostCommunication )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate &
      ( Communicator  =>  C % Manifold % Communicator, &
        nCB  =>  C % nCellsBrick, &
        nGL  =>  C % nGhostLayers, &
        nD   =>  C % nDimensions )
        
    select case ( nD )
    case ( 1 ) 
      return
    case ( 2 )
      DimensionMask = [ .true., .false., .false. ]
    case ( 3 )
      DimensionMask = [ .true., .true., .true. ]
    end select !-- nD

    !-- Allocate on First use
    
    if ( .not. allocated ( IncomingEdge ) &
         .and. .not. allocated ( OutgoingEdge ) ) then
    
      allocate ( IncomingEdge )
      allocate ( OutgoingEdge )
      
      !-- Post Receives

      call IncomingEdge % Initialize &
             ( Communicator, pack ( TagReceive, DimensionMask ), PH % Source, &
               PH % nChunksFrom  *  FSC % nFields )           
      call OutgoingEdge % Initialize &
             ( Communicator, pack ( TagSend, DimensionMask ), PH % Target, &
               PH % nChunksTo  *  FSC % nFields )
      
      if ( FSC % UseDeviceGhost ) then
        call IncomingEdge % AllocateDevice ( )
        call OutgoingEdge % AllocateDevice ( )
      end if
    
    end if  !-- allocated edges
      
    call T % Start ( )
    call IncomingEdge % Receive ( )
    call T % Stop ( )
    
    !-- Post Sends

    do kD = 3, 1, -1

      iD  =  mod ( kD, 3 ) + 1
      jD  =  mod ( iD, 3 ) + 1

      if ( iD > nD .or. jD > nD ) &
        cycle

      nSend ( iD )  =  nGL ( iD )
      nSend ( jD )  =  nGL ( jD )
      nSend ( kD )  =  nCB ( kD )

      !-- In setting oSend, note Copy command does not inherit lbound
      if ( TagSend ( kD )  ==  TAG_SEND_EDGE_RR ( kD ) ) then
        oSend         =  nGL
        oSend ( iD )  =  oSend ( iD ) + nCB ( iD ) - nGL ( iD )
        oSend ( jD )  =  oSend ( jD ) + nCB ( jD ) - nGL ( jD )
      else if ( TagSend ( kD )  ==  TAG_SEND_EDGE_LL ( kD ) ) then
        oSend         =  nGL
      else if ( TagSend ( kD )  ==  TAG_SEND_EDGE_RL ( kD ) ) then
        oSend         =  nGL
        oSend ( iD )  =  oSend ( iD ) + nCB ( iD ) - nGL ( iD )
      else if ( TagSend ( kD )  ==  TAG_SEND_EDGE_LR ( kD ) ) then
        oSend         =  nGL
        oSend ( jD )  =  oSend ( jD ) + nCB ( jD ) - nGL ( jD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'Field_CB__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartExchangeEdge', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      select case ( nD )
      case ( 2 )
        kM = 1
      case ( 3 )
        kM = kD
      end select !-- nD

      call LoadMessage &
             ( FSC, OutgoingEdge % Message ( kM ), nSend, oSend )
      
      call T % Start ( )
      call OutgoingEdge % Send ( kM )
      call T % Stop ( )

    end do  !-- kD

    !-- Cleanup

    end associate  !-- Communicator, etc.
    end select  !-- C

    nullify ( T )

  end subroutine StartExchangeEdge


  subroutine FinishExchangeEdge &
               ( FSC, IncomingEdge, OutgoingEdge, TagReceive )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    type ( MessageIncoming_1D_R_Form ), intent ( inout ) :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      OutgoingEdge
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive

    integer ( KDI ) :: &
      kM, &   !-- kMessage
      iD, jD  !-- iDimension, etc.
    integer ( KDI ), dimension ( 3 ) :: &
      oReceive, &
      nReceive
    logical ( KDL ) :: &
      AllFinished
    type ( TimerForm ), pointer :: &
      T 
      
    T => PROGRAM_HEADER % TimerPointer ( FSC % iTimerGhostCommunication )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate &
      ( nCB  =>  C % nCellsBrick, &
        nGL  =>  C % nGhostLayers, &
         nD  =>  C % nDimensions, &
        iaB  =>  C % iaBrick, &
         nB  =>  C % nBricks )

    if ( nD == 1 ) &
      return

    !-- Wait for Receives

    do 

      call T % Start ( )
      call IncomingEdge % Wait ( AllFinished, kM )
      call T % Stop ( )
      
      if ( AllFinished ) exit

      select case ( nD )
      case ( 2 )
        iD = 1
        jD = 2
      case ( 3 )
        iD = mod ( kM, 3 ) + 1
        jD = mod ( iD, 3 ) + 1
      end select !-- nD

      nReceive         =  nCB
      nReceive ( iD )  =  nGL ( iD )
      nReceive ( jD )  =  nGL ( jD )

      !-- In setting oReceive, note Copy command does not inherit lbound
      if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_LL ( kM ) ) then
        if ( iaB ( iD )  ==  1  .and.  iaB ( jD ) == 1  &
             .and..not. C % Periodic ( iD ) .and..not. C % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  -  nGL ( iD )
        oReceive ( jD )  =  oReceive ( jD )  -  nGL ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RR ( kM ) ) then
        if ( iaB ( iD )  ==  nB ( iD )  .and.  iaB ( jD )  ==  nB ( jD )  &
             .and..not. C % Periodic ( iD ) .and..not. C % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
        oReceive ( jD )  =  oReceive ( jD )  +  nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_LR ( kM ) ) then
        if ( iaB ( iD )  ==  1  .and.  iaB ( jD )  ==  nB ( jD )  &
             .and..not. C % Periodic ( iD ) .and..not. C % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  -  nGL ( iD )
        oReceive ( jD )  =  oReceive ( jD )  +  nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RL ( kM ) ) then
        if ( iaB ( iD )  ==  nB ( iD )  .and.  iaB ( jD )  ==  1  &
             .and. .not. C % Periodic ( iD ) &
             .and. .not. C % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
        oReceive ( jD )  =  oReceive ( jD )  -  nGL ( jD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_CB__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchangeEdge', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage &
             ( FSC, IncomingEdge % Message ( kM ), nReceive, oReceive )

    end do

    !-- Wait for Sends
    call T % Start ( )
    call OutgoingEdge % Wait ( )
    call T % Stop ( )

    !-- Cleanup

    end associate  !-- nCB, etc.
    end select  !-- C
    
    nullify ( T )

  end subroutine FinishExchangeEdge


  subroutine LoadMessage ( FSC, OutgoingMessage, nSend, oSend )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    type ( MessageOutgoing_R_Form ), intent ( in ) :: &
      OutgoingMessage
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nSend, &
      oSend

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F  !-- Field
    type ( TimerForm ), pointer :: &
      T

    T => PROGRAM_HEADER % TimerPointer ( FSC % iTimerGhostPackUnpack )
    call T % Start ( )
    
    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate ( FS  =>  FSC % FieldSet )

    oBuffer = 0
    do iS = 1, FS % nVariables
      iF = FS % iaSelected ( iS )
      call C % SetFieldPointer ( FS % Value ( :, iF ), F )
      call Copy ( F, nSend, oSend, oBuffer, OutgoingMessage % Value, &
                  UseDeviceOption = FSC % UseDeviceGhost )
      oBuffer = oBuffer + product ( nSend )
    end do !-- iS

    end associate !-- FS
    end select !-- C
    nullify ( F )

    call T % Stop ( )
    nullify ( T )

  end subroutine LoadMessage


  subroutine StoreMessage ( FSC, IncomingMessage, nReceive, oReceive )
               
    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC
    type ( MessageIncoming_R_Form ), intent ( in ) :: &
      IncomingMessage
    integer ( KDI ), dimension ( 3 ), intent ( in )  :: &
      nReceive, &
      oReceive

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F  !-- Field
    type ( TimerForm ), pointer :: &
      T
    
    T => PROGRAM_HEADER % TimerPointer ( FSC % iTimerGhostPackUnpack )
    call T % Start ( )
    
    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate ( FS  =>  FSC % FieldSet )

    oBuffer = 0
    do iS = 1, FS % nVariables          
      iF = FS % iaSelected ( iS )
      call C % SetFieldPointer ( FS % Value ( :, iF ), F )
      call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, F, &
                  UseDeviceOption = FSC % UseDeviceGhost )
      oBuffer = oBuffer + product ( nReceive )
    end do !-- iS
    
    end associate !-- FS
    end select !-- C
    nullify ( F )

    call T % Stop ( )    
    nullify ( T )

  end subroutine StoreMessage


end module FieldSet_CB__Form
