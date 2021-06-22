module GhostExchange_FSC__Form

  !-- GhostExchange_FieldSetChart_Form

  use Basics
  use Manifolds
  use Storage_FSC__Form

  implicit none
  private

  type, public :: GhostExchange_FSC_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer_C  = 0, &  !-- Communication
      iTimer_PU = 0     !-- PackUnpack
    logical ( KDL ) :: &
      DevicesCommunicate
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
      Exchange
    procedure, public, pass :: &
      StartExchange
    procedure, public, pass :: &
      FinishExchange
    final :: &
      Finalize
  end type GhostExchange_FSC_Form

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


  subroutine Initialize ( GE, DevicesCommunicateOption )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    logical ( KDL ), intent ( in ), optional :: &
      DevicesCommunicateOption

    GE % IGNORABILITY  =  CONSOLE % INFO_4

    GE % DevicesCommunicate  =  .false.
    if ( present ( DevicesCommunicateOption ) )  &
      GE % DevicesCommunicate  =  DevicesCommunicateOption  

  end subroutine Initialize


  subroutine Exchange  ( GE, SFSC, C, TimerLevelOption )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    class ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call GE % StartExchange ( C, SFSC, TimerLevelOption )
    call GE % FinishExchange ( SFSC, C )
   
  end subroutine Exchange


  subroutine StartExchange  ( GE, C, SFSC, TimerLevelOption )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    class ( Storage_FSC_Form ), intent ( in ) :: &
      SFSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDF ) :: &
      TimerName

    associate ( S  =>  SFSC % Storage )

    call Show ( 'Starting ghost exchange', GE % IGNORABILITY )
    call Show ( S % Name, 'FieldSet', GE % IGNORABILITY )

    associate &
      ( iT_GC   =>  GE % iTimer_C, &
        iT_GPU  =>  GE % iTimer_PU )
    if ( iT_GC == 0 ) then
      TimerName  =  'G_C_' // trim ( S % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GC, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GC, Level = 1 )
      end if
    end if
    if ( iT_GPU == 0 ) then
      TimerName  =  'G_PU_' // trim ( S % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GPU, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GPU, Level = 1 )
      end if
    end if
    end associate !-- iT_GC, etc.
      
    end associate !-- S

    select type ( C )
    class is ( Chart_GS_Form )

      associate ( S  =>  SFSC % Storage )

      !-- Start faces
      call StartFace_CGS &
             ( GE, GE % IncomingFace_L_R, GE % OutgoingFace_L_R, &
               C, C % PortalFace_L_R, S, TAG_RECEIVE_FACE_L, TAG_SEND_FACE_R )
      call StartFace_CGS &
             ( GE, GE % IncomingFace_R_L, GE % OutgoingFace_R_L, &
               C, C % PortalFace_R_L, S, TAG_RECEIVE_FACE_R, TAG_SEND_FACE_L )

      !-- Start edges
      call StartEdge_CGS &
             ( GE, GE % IncomingEdge_LL_RR, GE % OutgoingEdge_LL_RR, &
               C, C % PortalEdge_LL_RR, S, &
               TAG_RECEIVE_EDGE_LL, TAG_SEND_EDGE_RR )
      call StartEdge_CGS &
             ( GE, GE % IncomingEdge_RR_LL, GE % OutgoingEdge_RR_LL, &
               C, C % PortalEdge_RR_LL, S, &
               TAG_RECEIVE_EDGE_RR, TAG_SEND_EDGE_LL )
      call StartEdge_CGS &
             ( GE, GE % IncomingEdge_LR_RL, GE % OutgoingEdge_LR_RL, &
               C, C % PortalEdge_LR_RL, S, &
               TAG_RECEIVE_EDGE_LR, TAG_SEND_EDGE_RL )
      call StartEdge_CGS &
             ( GE, GE % IncomingEdge_RL_LR, GE % OutgoingEdge_RL_LR, &
               C, C % PortalEdge_RL_LR, S, &
               TAG_RECEIVE_EDGE_RL, TAG_SEND_EDGE_LR )

      end associate !-- S

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'GhostExchange_FSC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'StartExchange', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select  !-- C

  end subroutine StartExchange


  subroutine FinishExchange  ( GE, SFSC, C )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    class ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C

    associate ( S  =>  SFSC % Storage )
    call Show ( 'Finishing ghost exchange', GE % IGNORABILITY )
    call Show ( S % Name, 'FieldSet', GE % IGNORABILITY )
    end associate !-- S

    select type ( C )
    class is ( Chart_GS_Form )

      associate ( S  =>  SFSC % Storage )

      !-- Finish faces
      call FinishFace_CGS &
             ( GE, GE % IncomingFace_L_R, GE % OutgoingFace_L_R, S, C, &
               TAG_RECEIVE_FACE_L )
      call FinishFace_CGS &
             ( GE, GE % IncomingFace_R_L, GE % OutgoingFace_R_L, S, C, &
               TAG_RECEIVE_FACE_R )

      !-- Finish edges
      call FinishEdge_CGS &
             ( GE, GE % IncomingEdge_LL_RR, GE % OutgoingEdge_LL_RR, S, C, &
               TAG_RECEIVE_EDGE_LL )
      call FinishEdge_CGS &
             ( GE, GE % IncomingEdge_RR_LL, GE % OutgoingEdge_RR_LL, S, C, &
               TAG_RECEIVE_EDGE_RR )
      call FinishEdge_CGS &
             ( GE, GE % IncomingEdge_LR_RL, GE % OutgoingEdge_LR_RL, S, C, &
               TAG_RECEIVE_EDGE_LR )
      call FinishEdge_CGS &
             ( GE, GE % IncomingEdge_RL_LR, GE % OutgoingEdge_RL_LR, S, C, &
               TAG_RECEIVE_EDGE_RL )

      end associate !-- S

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'GhostExchange_FSC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'FinishExchange', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select  !-- C

  end subroutine FinishExchange


  impure elemental subroutine Finalize ( GE )

    type ( GhostExchange_FSC_Form ), intent ( inout ) :: &
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
               ( GE, IncomingFace, OutgoingFace, C, PH, S, TagReceive, TagSend )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
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

    T  =>  PROGRAM_HEADER % TimerPointer ( GE % iTimer_C )

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
    
      if ( GE % DevicesCommunicate ) then
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
        call Show ( 'FieldSet_CGS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartFace_CGS', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      call LoadMessage_CGS &
             ( GE, C, S, OutgoingFace % Message ( iD ), nSend, oSend )
      
      call T % Start ( )
      call OutgoingFace % Send ( iD )
      call T % Stop ( )

    end do !-- iD

    !-- Cleanup

    end associate  !-- Communicator, etc.

    nullify ( T )

  end subroutine StartFace_CGS


  subroutine FinishFace_CGS &
               ( GE, IncomingFace, OutgoingFace, S, C, TagReceive )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    type ( MessageIncoming_1D_R_Form ), intent ( inout ) :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      OutgoingFace
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
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
      
    T => PROGRAM_HEADER % TimerPointer ( GE % iTimer_C )

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
        call Show ( 'FieldSet_CGS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage_CGS &
             ( GE, S, C, IncomingFace % Message ( iD ), nReceive, oReceive )

    end do

    !-- Wait for Sends
    call T % Start ( )
    call OutgoingFace % Wait ( )
    call T % Stop ( )

    !-- Cleanup

    end associate !-- nCB etc.
    
    nullify ( T )

  end subroutine FinishFace_CGS


  subroutine StartEdge_CGS &
               ( GE, IncomingEdge, OutgoingEdge, C, PH, S, TagReceive, TagSend )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
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
      
    T => PROGRAM_HEADER % TimerPointer ( GE % iTimer_C )

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
      
      if ( GE % DevicesCommunicate ) then
        call IncomingEdge % AllocateDevice ( )
        call OutgoingEdge % AllocateDevice ( )
      end if
    
    end if  !-- allocated edges
      
    !-- Post Receives

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
             ( GE, C, S, OutgoingEdge % Message ( kM ), nSend, oSend )
      
      call T % Start ( )
      call OutgoingEdge % Send ( kM )
      call T % Stop ( )

    end do  !-- kD

    !-- Cleanup

    end associate  !-- Communicator, etc.

    nullify ( T )

  end subroutine StartEdge_CGS


  subroutine FinishEdge_CGS &
               ( GE, IncomingEdge, OutgoingEdge, S, C, TagReceive )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    type ( MessageIncoming_1D_R_Form ), intent ( inout ) :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ) :: &
      OutgoingEdge
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
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
      
    T => PROGRAM_HEADER % TimerPointer ( GE % iTimer_C )

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
        call Show ( 'FieldSet_CGS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishEdge_CGS', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage_CGS &
             ( GE, S, C, IncomingEdge % Message ( kM ), nReceive, oReceive )

    end do

    !-- Wait for Sends
    call T % Start ( )
    call OutgoingEdge % Wait ( )
    call T % Stop ( )

    !-- Cleanup

    end associate  !-- nCB, etc.
    
    nullify ( T )

  end subroutine FinishEdge_CGS


  subroutine LoadMessage_CGS ( GE, C, S, OutgoingMessage, nSend, oSend )

    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    class ( StorageForm ), intent ( in ) :: &
      S
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

    T => PROGRAM_HEADER % TimerPointer ( GE % iTimer_PU )
    call T % Start ( )
    
    oBuffer = 0
    do iS = 1, S % nVariables
      iF = S % iaSelected ( iS )
      call C % SetFieldPointer ( S % Value ( :, iF ), F )
      call Copy ( F, nSend, oSend, oBuffer, OutgoingMessage % Value, &
                  UseDeviceOption = GE % DevicesCommunicate )
      oBuffer = oBuffer + product ( nSend )
    end do !-- iS
    nullify ( F )

    call T % Stop ( )
    nullify ( T )

  end subroutine LoadMessage_CGS


  subroutine StoreMessage_CGS ( GE, S, C, IncomingMessage, nReceive, oReceive )
               
    class ( GhostExchange_FSC_Form ), intent ( inout ) :: &
      GE
    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
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
    
    T => PROGRAM_HEADER % TimerPointer ( GE % iTimer_PU )
    call T % Start ( )
    
    oBuffer = 0
    do iS = 1, S % nVariables          
      iF = S % iaSelected ( iS )
      call C % SetFieldPointer ( S % Value ( :, iF ), F )
      call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, F, &
                  UseDeviceOption = GE % DevicesCommunicate )
      oBuffer = oBuffer + product ( nReceive )
    end do !-- iS    
    nullify ( F )

    call T % Stop ( )    
    nullify ( T )

  end subroutine StoreMessage_CGS


end module GhostExchange_FSC__Form
