module FieldSet_GS__Form

  !-- FieldSet_GridStructured__Form

  use Basics
  use Manifolds
  use FieldSet_CH__Form

  implicit none
  private

  type, public, extends ( FieldSet_CH_Form ) :: FieldSet_GS_Form
    integer ( KDI ) :: &
      iTimerGhostCommunication = 0, &
      iTimerGhostPackUnpack    = 0, &
      iTimerUpdateDevice       = 0, &
      iTimerUpdateHost         = 0
    logical ( KDL ) :: &
      DeviceMemory, &
      PinnedMemory, &
      DevicesCommunicate
    class ( StorageForm ), allocatable :: &
      FieldSet
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
    procedure, private, pass :: &
      InitializeAllocate_GS
    generic, public :: &
      Initialize => InitializeAllocate_GS
    procedure, public, pass :: &
      Clone
    procedure, public, pass :: &
      ExchangeGhostData
    procedure, public, pass :: &
      StartGhostExchange
    procedure, public, pass :: &
      FinishGhostExchange
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_FS
    procedure, public, pass :: &
      UpdateHost => UpdateHost_FS
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateFieldSet
    procedure, private, pass :: &
      CloneFieldSet
  end type FieldSet_GS_Form


    private :: &
      StartExchangeFace, &
      FinishExchangeFace, &
      StartExchangeEdge, &
      FinishExchangeEdge

      private :: &
        LoadMessage, &
        StoreMessage

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


  subroutine InitializeAllocate_GS &
               ( FSG, G, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
    class ( Grid_S_Form ), intent ( inout ), target :: &
      G
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    if ( FSG % Type == '' ) &
      FSG % Type  =  'a FieldSet_GS' 

    call FSG % Initialize_H &
           ( G, FieldOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption, nFieldsOption )

    FSG % DeviceMemory  =  .false.
    if ( present ( DeviceMemoryOption ) ) &
      FSG % DeviceMemory  =  DeviceMemoryOption
    
    FSG % PinnedMemory  =  .false.
    if ( present ( PinnedMemoryOption ) ) &
      FSG % PinnedMemory  =  PinnedMemoryOption
    
    FSG % DevicesCommunicate  =  .false.
    if ( present ( DevicesCommunicateOption ) )  &
      FSG % DevicesCommunicate  =  DevicesCommunicateOption  

    call FSG % AllocateFieldSet ( )

  end subroutine InitializeAllocate_GS


  subroutine Clone ( FSC_T, FSC_S, NameOption, iaSelectedOption )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSC_T  !-- FSC_Target
    class ( FieldSet_CH_Form ), intent ( in ), target :: &
      FSC_S  !-- FSC_Source
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    if ( FSC_T % Type == '' ) &
      FSC_T % Type  =  'a FieldSet_GS' 

    call FSC_T % FieldSet_CH_Form % Clone &
           ( FSC_S, NameOption, iaSelectedOption )

    select type ( FSC_S )
    class is ( FieldSet_GS_Form )

    FSC_T % DeviceMemory        =  FSC_S % DeviceMemory    
    FSC_T % PinnedMemory        =  FSC_S % PinnedMemory
    FSC_T % DevicesCommunicate  =  FSC_S % DevicesCommunicate

    call FSC_T % CloneFieldSet ( FSC_S )

    end select !-- FSC_S

  end subroutine Clone


  subroutine ExchangeGhostData  ( FSG, TimerLevelOption )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call FSG % StartGhostExchange ( TimerLevelOption )
    call FSG % FinishGhostExchange ( )
   
  end subroutine ExchangeGhostData


  subroutine StartGhostExchange  ( FSG, TimerLevelOption )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDF ) :: &
      TimerName

    call Show ( 'Starting ghost exchange', FSG % IGNORABILITY + 2 )
    call Show ( FSG % Name, 'FieldSet', FSG % IGNORABILITY + 2 )

    associate &
      ( iT_GC   =>  FSG % iTimerGhostCommunication, &
        iT_GPU  =>  FSG % iTimerGhostPackUnpack )
    if ( iT_GC == 0 ) then
      TimerName  =  'GhostCommunication ' // trim ( FSG % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GC, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GC, Level = 1 )
      end if
    end if
    if ( iT_GPU == 0 ) then
      TimerName  =  'GhostPackUnpack ' // trim ( FSG % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GPU, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GPU, Level = 1 )
      end if
    end if
    end associate !-- iT_GC, etc.
      
    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    !-- Start faces
    call StartExchangeFace &
           ( FSG, FSG % IncomingFace_L_R, FSG % OutgoingFace_L_R, &
             G % PortalFace_L_R, TAG_RECEIVE_FACE_L, TAG_SEND_FACE_R )
    call StartExchangeFace &
           ( FSG, FSG % IncomingFace_R_L, FSG % OutgoingFace_R_L, &
             G % PortalFace_R_L, TAG_RECEIVE_FACE_R, TAG_SEND_FACE_L )

    !-- Start edges
    call StartExchangeEdge &
           ( FSG, FSG % IncomingEdge_LL_RR, FSG % OutgoingEdge_LL_RR, &
             G % PortalEdge_LL_RR, TAG_RECEIVE_EDGE_LL, TAG_SEND_EDGE_RR )
    call StartExchangeEdge &
           ( FSG, FSG % IncomingEdge_RR_LL, FSG % OutgoingEdge_RR_LL, &
             G % PortalEdge_RR_LL, TAG_RECEIVE_EDGE_RR, TAG_SEND_EDGE_LL )
    call StartExchangeEdge &
           ( FSG, FSG % IncomingEdge_LR_RL, FSG % OutgoingEdge_LR_RL, &
             G % PortalEdge_LR_RL, TAG_RECEIVE_EDGE_LR, TAG_SEND_EDGE_RL )
    call StartExchangeEdge &
           ( FSG, FSG % IncomingEdge_RL_LR, FSG % OutgoingEdge_RL_LR, &
             G % PortalEdge_RL_LR, TAG_RECEIVE_EDGE_RL, TAG_SEND_EDGE_LR )

    end select  !-- C

  end subroutine StartGhostExchange


  subroutine FinishGhostExchange  ( FSG )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG

    call Show ( 'Finishing ghost exchange', FSG % IGNORABILITY + 2 )
    call Show ( FSG % Name, 'FieldSet', FSG % IGNORABILITY + 2 )

    !-- Finish faces
    call FinishExchangeFace &
           ( FSG, FSG % IncomingFace_L_R, FSG % OutgoingFace_L_R, &
             TAG_RECEIVE_FACE_L )
    call FinishExchangeFace &
           ( FSG, FSG % IncomingFace_R_L, FSG % OutgoingFace_R_L, &
             TAG_RECEIVE_FACE_R )

    !-- Finish edges
    call FinishExchangeEdge &
           ( FSG, FSG % IncomingEdge_LL_RR, FSG % OutgoingEdge_LL_RR, &
             TAG_RECEIVE_EDGE_LL )
    call FinishExchangeEdge &
           ( FSG, FSG % IncomingEdge_RR_LL, FSG % OutgoingEdge_RR_LL, &
             TAG_RECEIVE_EDGE_RR )
    call FinishExchangeEdge &
           ( FSG, FSG % IncomingEdge_LR_RL, FSG % OutgoingEdge_LR_RL, &
             TAG_RECEIVE_EDGE_LR )
    call FinishExchangeEdge &
           ( FSG, FSG % IncomingEdge_RL_LR, FSG % OutgoingEdge_RL_LR, &
             TAG_RECEIVE_EDGE_RL )
    
  end subroutine FinishGhostExchange


  subroutine UpdateDevice_FS ( FSG, TimerLevelOption )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  FSG % iTimerUpdateDevice )
    if ( iT == 0 ) then
      TimerName  =  'UpdateDevice ' // trim ( FSG % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( FSG % iTimerUpdateDevice )

    call T % Start ( )
    call FSG % FieldSet % UpdateDevice ( )
    call T % Stop ( )

  end subroutine UpdateDevice_FS


  subroutine UpdateHost_FS ( FSG, TimerLevelOption )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  FSG % iTimerUpdateHost )
    if ( iT == 0 ) then
      TimerName  =  'UpdateHost ' // trim ( FSG % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( FSG % iTimerUpdateHost )

    call T % Start ( )
    call FSG % FieldSet % UpdateHost ( )
    call T % Stop ( )

  end subroutine UpdateHost_FS


  impure elemental subroutine Finalize ( FSG )

    type ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG

    if ( allocated ( FSG % OutgoingEdge_RL_LR ) ) &
      deallocate ( FSG % OutgoingEdge_RL_LR )
    if ( allocated ( FSG % OutgoingEdge_LR_RL ) ) &
      deallocate ( FSG % OutgoingEdge_LR_RL )
    if ( allocated ( FSG % OutgoingEdge_RR_LL ) ) &
      deallocate ( FSG % OutgoingEdge_RR_LL )
    if ( allocated ( FSG % OutgoingEdge_LL_RR ) ) &
      deallocate ( FSG % OutgoingEdge_LL_RR )
    if ( allocated ( FSG % OutgoingFace_R_L ) ) &
      deallocate ( FSG % OutgoingFace_R_L )
    if ( allocated ( FSG % OutgoingFace_L_R ) ) &
      deallocate ( FSG % OutgoingFace_L_R )

    if ( allocated ( FSG % IncomingEdge_RL_LR ) ) &
      deallocate ( FSG % IncomingEdge_RL_LR )
    if ( allocated ( FSG % IncomingEdge_LR_RL ) ) &
      deallocate ( FSG % IncomingEdge_LR_RL )
    if ( allocated ( FSG % IncomingEdge_RR_LL ) ) &
      deallocate ( FSG % IncomingEdge_RR_LL )
    if ( allocated ( FSG % IncomingEdge_LL_RR ) ) &
      deallocate ( FSG % IncomingEdge_LL_RR )
    if ( allocated ( FSG % IncomingFace_R_L ) ) &
      deallocate ( FSG % IncomingFace_R_L )
    if ( allocated ( FSG % IncomingFace_L_R ) ) &
      deallocate ( FSG % IncomingFace_L_R )

    if ( allocated ( FSG % FieldSet ) ) &
      deallocate ( FSG % FieldSet )

  end subroutine Finalize


  subroutine AllocateFieldSet ( FSG )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG

    select type ( G => FSG % Chart )
    class is ( Grid_S_Form )

    if ( .not. allocated ( FSG % FieldSet ) ) then
      allocate ( FSG % FieldSet )
      associate ( FS  =>  FSG % FieldSet )
      call FS % Initialize &
             ( [ G % nCellsLocal, FSG % nFields ], &
               VariableOption = FSG % Field, &
               VectorOption = FSG % Vector, &
               NameOption = FSG % Name, &
               ClearOption = .true., &
               PinnedOption = FSG % PinnedMemory, &
               UnitOption = FSG % Unit, &
               VectorIndicesOption = FSG % VectorIndices )
      if ( FSG % DeviceMemory ) &
        call FS % AllocateDevice ( )
      end associate !-- FS
    end if

    end select !-- G

  end subroutine AllocateFieldSet


  subroutine CloneFieldSet ( FSG_T, FSG_S )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG_T
    class ( FieldSet_GS_Form ), intent ( in ) :: &
      FSG_S

    if ( .not. allocated ( FSG_T % FieldSet ) ) then
      allocate ( FSG_T % FieldSet )
      associate ( FS  =>  FSG_T % FieldSet )
      call FS % Initialize &
             ( FSG_S % FieldSet, &
               VectorOption = FSG_T % Vector, &
               NameOption = FSG_T % Name, &
               VectorIndicesOption = FSG_T % VectorIndices, &
               iaSelectedOption = FSG_T % iaSelected )
      end associate !-- FS
    end if

  end subroutine CloneFieldSet


  subroutine StartExchangeFace &
               ( FSG, IncomingFace, OutgoingFace, PH, TagReceive, TagSend )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
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

    T  =>  PROGRAM_HEADER % TimerPointer ( FSG % iTimerGhostCommunication )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate &
      ( Communicator  =>  G % Communicator, &
        nCB  =>  G % nCellsBrick, &
        nGL  =>  G % nGhostLayers, &
        nD   =>  G % nDimensions )

    !-- Allocate on first use

    if ( .not. allocated ( IncomingFace ) &
         .and. .not. allocated ( OutgoingFace ) ) then
    
      allocate ( IncomingFace )
      allocate ( OutgoingFace )

      call IncomingFace % Initialize &
             ( Communicator, TagReceive ( : nD ), PH % Source, &
               PH % nChunksFrom  *  FSG % nFields )
      call OutgoingFace % Initialize &
             ( Communicator, TagSend ( : nD ), PH % Target, &
               PH % nChunksTo  *  FSG % nFields )
    
      if ( FSG % DevicesCommunicate ) then
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
        call Show ( 'FieldSet_GS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartExchangeFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      call LoadMessage &
             ( FSG, OutgoingFace % Message ( iD ), nSend, oSend )
      
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
               ( FSG, IncomingFace, OutgoingFace, TagReceive )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
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
      
    T => PROGRAM_HEADER % TimerPointer ( FSG % iTimerGhostCommunication )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate &
      ( nCB => G % nCellsBrick, &
        nGL => G % nGhostLayers, &
        iaB => G % iaBrick, &
         nB => G % nBricks )

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
        if ( iaB ( iD )  ==  1 .and. .not. G % Periodic ( iD ) ) &
          cycle
        oReceive        =  nGL
        oReceive ( iD ) =  oReceive ( iD )  -  nGL ( iD )
      else if ( TagReceive ( iD )  ==  TAG_RECEIVE_FACE_R ( iD ) ) then
        if ( iaB ( iD )  ==  nB ( iD ) .and. .not. G % Periodic ( iD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_GS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchangeFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage &
             ( FSG, IncomingFace % Message ( iD ), nReceive, oReceive )

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
               ( FSG, IncomingEdge, OutgoingEdge, PH, TagReceive, TagSend )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
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
      
    T => PROGRAM_HEADER % TimerPointer ( FSG % iTimerGhostCommunication )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate &
      ( Communicator  =>  G % Communicator, &
        nCB  =>  G % nCellsBrick, &
        nGL  =>  G % nGhostLayers, &
        nD   =>  G % nDimensions )
        
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
               PH % nChunksFrom  *  FSG % nFields )           
      call OutgoingEdge % Initialize &
             ( Communicator, pack ( TagSend, DimensionMask ), PH % Target, &
               PH % nChunksTo  *  FSG % nFields )
      
      if ( FSG % DevicesCommunicate ) then
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
        call Show ( 'Field_GS__Form', 'module', CONSOLE % ERROR )
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
             ( FSG, OutgoingEdge % Message ( kM ), nSend, oSend )
      
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
               ( FSG, IncomingEdge, OutgoingEdge, TagReceive )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
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
      
    T => PROGRAM_HEADER % TimerPointer ( FSG % iTimerGhostCommunication )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate &
      ( nCB  =>  G % nCellsBrick, &
        nGL  =>  G % nGhostLayers, &
         nD  =>  G % nDimensions, &
        iaB  =>  G % iaBrick, &
         nB  =>  G % nBricks )

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
             .and..not. G % Periodic ( iD ) .and..not. G % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  -  nGL ( iD )
        oReceive ( jD )  =  oReceive ( jD )  -  nGL ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RR ( kM ) ) then
        if ( iaB ( iD )  ==  nB ( iD )  .and.  iaB ( jD )  ==  nB ( jD )  &
             .and..not. G % Periodic ( iD ) .and..not. G % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
        oReceive ( jD )  =  oReceive ( jD )  +  nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_LR ( kM ) ) then
        if ( iaB ( iD )  ==  1  .and.  iaB ( jD )  ==  nB ( jD )  &
             .and..not. G % Periodic ( iD ) .and..not. G % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  -  nGL ( iD )
        oReceive ( jD )  =  oReceive ( jD )  +  nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RL ( kM ) ) then
        if ( iaB ( iD )  ==  nB ( iD )  .and.  iaB ( jD )  ==  1  &
             .and. .not. G % Periodic ( iD ) &
             .and. .not. G % Periodic ( jD ) ) &
          cycle
        oReceive         =  nGL
        oReceive ( iD )  =  oReceive ( iD )  +  nCB ( iD )
        oReceive ( jD )  =  oReceive ( jD )  -  nGL ( jD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'FieldSet_GS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchangeEdge', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage &
             ( FSG, IncomingEdge % Message ( kM ), nReceive, oReceive )

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


  subroutine LoadMessage ( FSG, OutgoingMessage, nSend, oSend )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
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

    T => PROGRAM_HEADER % TimerPointer ( FSG % iTimerGhostPackUnpack )
    call T % Start ( )
    
    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate ( FS  =>  FSG % FieldSet )

    oBuffer = 0
    do iS = 1, FS % nVariables
      iF = FS % iaSelected ( iS )
      call G % SetFieldPointer ( FS % Value ( :, iF ), F )
      call Copy ( F, nSend, oSend, oBuffer, OutgoingMessage % Value, &
                  UseDeviceOption = FSG % DevicesCommunicate )
      oBuffer = oBuffer + product ( nSend )
    end do !-- iS

    end associate !-- FS
    end select !-- C
    nullify ( F )

    call T % Stop ( )
    nullify ( T )

  end subroutine LoadMessage


  subroutine StoreMessage ( FSG, IncomingMessage, nReceive, oReceive )
               
    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
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
    
    T => PROGRAM_HEADER % TimerPointer ( FSG % iTimerGhostPackUnpack )
    call T % Start ( )
    
    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate ( FS  =>  FSG % FieldSet )

    oBuffer = 0
    do iS = 1, FS % nVariables          
      iF = FS % iaSelected ( iS )
      call G % SetFieldPointer ( FS % Value ( :, iF ), F )
      call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, F, &
                  UseDeviceOption = FSG % DevicesCommunicate )
      oBuffer = oBuffer + product ( nReceive )
    end do !-- iS
    
    end associate !-- FS
    end select !-- C
    nullify ( F )

    call T % Stop ( )    
    nullify ( T )

  end subroutine StoreMessage


end module FieldSet_GS__Form
