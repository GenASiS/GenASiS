module GhostExchange_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: GhostExchangeForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
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
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      StartExchange
    procedure, public, pass :: &
      FinishExchange
    final :: &
      Finalize
  end type GhostExchangeForm

    private :: &
      StartFace_CGS, &
      FinishFace_CGS, &
      StartEdge_CGS, &
      FinishEdge_CGS

      private :: &
        LoadMessage_CGS, &
        StoreMessage_CGS

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


  subroutine Initialize ( GE )

    class ( GhostExchangeForm ), intent ( inout ) :: &
      GE

    GE % IGNORABILITY  =  CONSOLE % INFO_4

  end subroutine Initialize


  subroutine StartExchange ( GE, C, S, DevicesCommunicate )

    class ( GhostExchangeForm ), intent ( inout ) :: &
      GE
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    class ( StorageForm ), intent ( in ) :: &
      S
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate

    call Show ( 'Starting ghost exchange', GE % IGNORABILITY )
    call Show ( S % Name, 'FieldSet', GE % IGNORABILITY )

    select type ( C )
    class is ( Chart_GS_Form )

      !-- Start faces
      call StartFace_CGS &
             ( GE % IncomingFace_L_R, GE % OutgoingFace_L_R, &
               C, C % PortalFace_L_R, S, DevicesCommunicate, &
               TAG_RECEIVE_FACE_L, TAG_SEND_FACE_R )
      call StartFace_CGS &
             ( GE % IncomingFace_R_L, GE % OutgoingFace_R_L, &
               C, C % PortalFace_R_L, S, DevicesCommunicate, &
               TAG_RECEIVE_FACE_R, TAG_SEND_FACE_L )

      !-- Start edges
      call StartEdge_CGS &
             ( GE % IncomingEdge_LL_RR, GE % OutgoingEdge_LL_RR, &
               C, C % PortalEdge_LL_RR, S, DevicesCommunicate, &
               TAG_RECEIVE_EDGE_LL, TAG_SEND_EDGE_RR )
      call StartEdge_CGS &
             ( GE % IncomingEdge_RR_LL, GE % OutgoingEdge_RR_LL, &
               C, C % PortalEdge_RR_LL, S, DevicesCommunicate, &
               TAG_RECEIVE_EDGE_RR, TAG_SEND_EDGE_LL )
      call StartEdge_CGS &
             ( GE % IncomingEdge_LR_RL, GE % OutgoingEdge_LR_RL, &
               C, C % PortalEdge_LR_RL, S, DevicesCommunicate, &
               TAG_RECEIVE_EDGE_LR, TAG_SEND_EDGE_RL )
      call StartEdge_CGS &
             ( GE % IncomingEdge_RL_LR, GE % OutgoingEdge_RL_LR, &
               C, C % PortalEdge_RL_LR, S, DevicesCommunicate, &
               TAG_RECEIVE_EDGE_RL, TAG_SEND_EDGE_LR )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'GhostExchangeForm', 'module', CONSOLE % ERROR )
      call Show ( 'StartExchange', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select  !-- C

  end subroutine StartExchange


  subroutine FinishExchange ( GE, S, C, Periodic, DevicesCommunicate )

    class ( GhostExchangeForm ), intent ( inout ) :: &
      GE
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      Periodic
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate

    call Show ( 'Finishing ghost exchange', GE % IGNORABILITY )
    call Show ( S % Name, 'FieldSet', GE % IGNORABILITY )

    select type ( C )
    class is ( Chart_GS_Form )

      !-- Finish faces
      call FinishFace_CGS &
             ( GE % IncomingFace_L_R, GE % OutgoingFace_L_R, S, C, &
               Periodic, DevicesCommunicate, TAG_RECEIVE_FACE_L )
      call FinishFace_CGS &
             ( GE % IncomingFace_R_L, GE % OutgoingFace_R_L, S, C, &
               Periodic, DevicesCommunicate, TAG_RECEIVE_FACE_R )

      !-- Finish edges
      call FinishEdge_CGS &
             ( GE % IncomingEdge_LL_RR, GE % OutgoingEdge_LL_RR, S, C, &
               Periodic, DevicesCommunicate, TAG_RECEIVE_EDGE_LL )
      call FinishEdge_CGS &
             ( GE % IncomingEdge_RR_LL, GE % OutgoingEdge_RR_LL, S, C, &
               Periodic, DevicesCommunicate, TAG_RECEIVE_EDGE_RR )
      call FinishEdge_CGS &
             ( GE % IncomingEdge_LR_RL, GE % OutgoingEdge_LR_RL, S, C, &
               Periodic, DevicesCommunicate, TAG_RECEIVE_EDGE_LR )
      call FinishEdge_CGS &
             ( GE % IncomingEdge_RL_LR, GE % OutgoingEdge_RL_LR, S, C, &
               Periodic, DevicesCommunicate, TAG_RECEIVE_EDGE_RL )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'GhostExchangeForm', 'module', CONSOLE % ERROR )
      call Show ( 'FinishExchange', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select  !-- C

  end subroutine FinishExchange


  impure elemental subroutine Finalize ( GE )

    type ( GhostExchangeForm ), intent ( inout ) :: &
      GE

    if ( allocated ( GE % OutgoingEdge_RL_LR ) ) &
      deallocate ( GE % OutgoingEdge_RL_LR )
    if ( allocated ( GE % OutgoingEdge_LR_RL ) ) &
      deallocate ( GE % OutgoingEdge_LR_RL )
    if ( allocated ( GE % OutgoingEdge_RR_LL ) ) &
      deallocate ( GE % OutgoingEdge_RR_LL )
    if ( allocated ( GE % OutgoingEdge_LL_RR ) ) &
      deallocate ( GE % OutgoingEdge_LL_RR )
    if ( allocated ( GE % OutgoingFace_R_L ) ) &
      deallocate ( GE % OutgoingFace_R_L )
    if ( allocated ( GE % OutgoingFace_L_R ) ) &
      deallocate ( GE % OutgoingFace_L_R )

    if ( allocated ( GE % IncomingEdge_RL_LR ) ) &
      deallocate ( GE % IncomingEdge_RL_LR )
    if ( allocated ( GE % IncomingEdge_LR_RL ) ) &
      deallocate ( GE % IncomingEdge_LR_RL )
    if ( allocated ( GE % IncomingEdge_RR_LL ) ) &
      deallocate ( GE % IncomingEdge_RR_LL )
    if ( allocated ( GE % IncomingEdge_LL_RR ) ) &
      deallocate ( GE % IncomingEdge_LL_RR )
    if ( allocated ( GE % IncomingFace_R_L ) ) &
      deallocate ( GE % IncomingFace_R_L )
    if ( allocated ( GE % IncomingFace_L_R ) ) &
      deallocate ( GE % IncomingFace_L_R )

  end subroutine Finalize


  subroutine StartFace_CGS &
               ( IncomingFace, OutgoingFace, C, PH, S, DevicesCommunicate, &
                 TagReceive, TagSend )

    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingFace
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    class ( StorageForm ), intent ( in ) :: &
      S
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend
    
    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend

    associate &
      ( Communicator  =>  C % Communicator, &
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
               PH % nChunksFrom  *  S % nVariables )
      call OutgoingFace % Initialize &
             ( Communicator, TagSend ( : nD ), PH % Target, &
               PH % nChunksTo  *  S % nVariables )
    
      if ( DevicesCommunicate ) then
        call IncomingFace % AllocateDevice ( )
        call OutgoingFace % AllocateDevice ( )
      end if 
    
    end if  !-- allocated faces
    
    !-- Post Receives

    call IncomingFace % Receive ( )

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
        call Show ( 'FieldSet_CGS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartFace_CGS', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      call LoadMessage_CGS &
             ( C, S, OutgoingFace % Message ( iD ), DevicesCommunicate, &
               nSend, oSend )
      
      call OutgoingFace % Send ( iD )

    end do !-- iD

    !-- Cleanup

    end associate  !-- Communicator, etc.

  end subroutine StartFace_CGS


  subroutine FinishFace_CGS &
               ( IncomingFace, OutgoingFace, S, C, Periodic, &
                 DevicesCommunicate, TagReceive )

    type ( MessageIncoming_1D_R_Form ), intent ( inout ) :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      OutgoingFace
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      Periodic
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive

    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oReceive, &
      nReceive
    logical ( KDL ) :: &
      AllFinished

    associate &
      ( nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        iaB => C % iaBrick, &
         nB => C % nBricks )

    !-- Wait for Receives

    do 

      call IncomingFace % Wait ( AllFinished, iD )
      
      if ( AllFinished ) exit

      nReceive        = nCB
      nReceive ( iD ) = nGL ( iD )

      !-- In setting oReceive, note Copy command does not inherit lbound
      if ( TagReceive ( iD )  ==  TAG_RECEIVE_FACE_L ( iD ) ) then
        if ( iaB ( iD )  ==  1 .and. .not. Periodic ( iD ) ) &
          cycle
        oReceive        =  nGL
        oReceive ( iD ) =  oReceive ( iD )  -  nGL ( iD )
      else if ( TagReceive ( iD )  ==  TAG_RECEIVE_FACE_R ( iD ) ) then
        if ( iaB ( iD )  ==  nB ( iD ) .and. .not. Periodic ( iD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_CGS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage_CGS &
             ( S, C, IncomingFace % Message ( iD ), DevicesCommunicate, &
               nReceive, oReceive )

    end do

    !-- Wait for Sends
    call OutgoingFace % Wait ( )

    !-- Cleanup

    end associate !-- nCB etc.
    
  end subroutine FinishFace_CGS


  subroutine StartEdge_CGS &
               ( IncomingEdge, OutgoingEdge, C, PH, S, DevicesCommunicate, &
                 TagReceive, TagSend )

    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingEdge
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    class ( StorageForm ), intent ( in ) :: &
      S
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
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

    associate &
      ( Communicator  =>  C % Communicator, &
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
      
      call IncomingEdge % Initialize &
             ( Communicator, pack ( TagReceive, DimensionMask ), PH % Source, &
               PH % nChunksFrom  *  S % nVariables )           
      call OutgoingEdge % Initialize &
             ( Communicator, pack ( TagSend, DimensionMask ), PH % Target, &
               PH % nChunksTo  *  S % nVariables )
      
      if ( DevicesCommunicate ) then
        call IncomingEdge % AllocateDevice ( )
        call OutgoingEdge % AllocateDevice ( )
      end if
    
    end if  !-- allocated edges
      
    !-- Post Receives

    call IncomingEdge % Receive ( )
    
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
        call Show ( 'Field_GS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartEdge_CGS', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      select case ( nD )
      case ( 2 )
        kM = 1
      case ( 3 )
        kM = kD
      end select !-- nD

      call LoadMessage_CGS &
             ( C, S, OutgoingEdge % Message ( kM ), DevicesCommunicate, &
               nSend, oSend )
      
      call OutgoingEdge % Send ( kM )

    end do  !-- kD

    !-- Cleanup

    end associate  !-- Communicator, etc.

  end subroutine StartEdge_CGS


  subroutine FinishEdge_CGS &
               ( IncomingEdge, OutgoingEdge, S, C, Periodic, &
                 DevicesCommunicate, TagReceive )

    type ( MessageIncoming_1D_R_Form ), intent ( inout ) :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      OutgoingEdge
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      Periodic
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
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

      call IncomingEdge % Wait ( AllFinished, kM )
      
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
             .and..not. Periodic ( iD ) .and..not. Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  -  nGL ( iD )
        oReceive ( jD )  =  oReceive ( jD )  -  nGL ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RR ( kM ) ) then
        if ( iaB ( iD )  ==  nB ( iD )  .and.  iaB ( jD )  ==  nB ( jD )  &
             .and..not. Periodic ( iD ) .and..not. Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
        oReceive ( jD )  =  oReceive ( jD )  +  nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_LR ( kM ) ) then
        if ( iaB ( iD )  ==  1  .and.  iaB ( jD )  ==  nB ( jD )  &
             .and..not. Periodic ( iD ) .and..not. Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  -  nGL ( iD )
        oReceive ( jD )  =  oReceive ( jD )  +  nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RL ( kM ) ) then
        if ( iaB ( iD )  ==  nB ( iD )  .and.  iaB ( jD )  ==  1  &
             .and. .not. Periodic ( iD ) &
             .and. .not. Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
        oReceive ( jD )  =  oReceive ( jD )  -  nGL ( jD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_CGS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishEdge_CGS', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage_CGS &
             ( S, C, IncomingEdge % Message ( kM ), DevicesCommunicate, &
               nReceive, oReceive )

    end do

    !-- Wait for Sends
    call OutgoingEdge % Wait ( )

    !-- Cleanup

    end associate  !-- nCB, etc.
    
  end subroutine FinishEdge_CGS


  subroutine LoadMessage_CGS &
               ( C, S, OutgoingMessage, DevicesCommunicate, nSend, oSend )

    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    class ( StorageForm ), intent ( in ) :: &
      S
    type ( MessageOutgoing_R_Form ), intent ( in ) :: &
      OutgoingMessage
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nSend, &
      oSend

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F  !-- Field

    oBuffer = 0
    do iS = 1, S % nVariables
      iF = S % iaSelected ( iS )
      call C % SetFieldPointer ( S % Value ( :, iF ), F )
      call Copy ( F, nSend, oSend, oBuffer, OutgoingMessage % Value, &
                  UseDeviceOption = DevicesCommunicate )
      oBuffer = oBuffer + product ( nSend )
    end do !-- iS
    nullify ( F )

  end subroutine LoadMessage_CGS


  subroutine StoreMessage_CGS &
               ( S, C, IncomingMessage, DevicesCommunicate, &
                 nReceive, oReceive )
               
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( MessageIncoming_R_Form ), intent ( in ) :: &
      IncomingMessage
    logical ( KDL ), intent ( in ) :: &
      DevicesCommunicate
    integer ( KDI ), dimension ( 3 ), intent ( in )  :: &
      nReceive, &
      oReceive

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F  !-- Field
    
    oBuffer = 0
    do iS = 1, S % nVariables          
      iF = S % iaSelected ( iS )
      call C % SetFieldPointer ( S % Value ( :, iF ), F )
      call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, F, &
                  UseDeviceOption = DevicesCommunicate )
      oBuffer = oBuffer + product ( nReceive )
    end do !-- iS    
    nullify ( F )

  end subroutine StoreMessage_CGS


end module GhostExchange_Form
