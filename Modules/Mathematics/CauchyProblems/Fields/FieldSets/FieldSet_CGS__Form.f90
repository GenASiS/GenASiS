module FieldSet_CGS__Form

  !-- FieldSet_GridStructured__Form

  use Basics
  use Manifolds
  use FieldSet_CH__Form

  implicit none
  private

  type, public, extends ( FieldSet_CH_Form ) :: FieldSet_CGS_Form
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
      Storage
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
      AllocateStorage
    procedure, private, pass :: &
      CloneStorage
  end type FieldSet_CGS_Form


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


  subroutine Initialize &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_GS_Form ), intent ( inout ), target :: &
      C
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

    if ( FSC % Type == '' ) &
      FSC % Type  =  'a FieldSet_CGS' 

    call FSC % Initialize_H &
           ( C, FieldOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption, nFieldsOption )

    FSC % DeviceMemory  =  .false.
    if ( present ( DeviceMemoryOption ) ) &
      FSC % DeviceMemory  =  DeviceMemoryOption
    
    FSC % PinnedMemory  =  .false.
    if ( present ( PinnedMemoryOption ) ) &
      FSC % PinnedMemory  =  PinnedMemoryOption
    
    FSC % DevicesCommunicate  =  .false.
    if ( present ( DevicesCommunicateOption ) )  &
      FSC % DevicesCommunicate  =  DevicesCommunicateOption  

    call FSC % AllocateStorage ( )

  end subroutine Initialize


  subroutine Clone &
               ( FSC_T, FSC_S, NameOption, iaSelectedOption, &
                 IgnorabilityOption )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC_T  !-- FSC_Target
    class ( FieldSet_CH_Form ), intent ( in ), target :: &
      FSC_S  !-- FSC_Source
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FSC_T % Type == '' ) &
      FSC_T % Type  =  'a FieldSet_CGS' 

    call FSC_T % FieldSet_CH_Form % Clone &
           ( FSC_S, NameOption, iaSelectedOption, IgnorabilityOption )

    select type ( FSC_S )
    class is ( FieldSet_CGS_Form )

    FSC_T % DeviceMemory        =  FSC_S % DeviceMemory    
    FSC_T % PinnedMemory        =  FSC_S % PinnedMemory
    FSC_T % DevicesCommunicate  =  FSC_S % DevicesCommunicate

    call FSC_T % CloneStorage ( FSC_S )

    end select !-- FSC_S

  end subroutine Clone


  subroutine ExchangeGhostData  ( FSC, TimerLevelOption )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call FSC % StartGhostExchange ( TimerLevelOption )
    call FSC % FinishGhostExchange ( )
   
  end subroutine ExchangeGhostData


  subroutine StartGhostExchange  ( FSC, TimerLevelOption )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDF ) :: &
      TimerName

    call Show ( 'Starting ghost exchange', FSC % IGNORABILITY + 2 )
    call Show ( FSC % Name, 'FieldSet', FSC % IGNORABILITY + 2 )

    associate &
      ( iT_GC   =>  FSC % iTimerGhostCommunication, &
        iT_GPU  =>  FSC % iTimerGhostPackUnpack )
    if ( iT_GC == 0 ) then
      TimerName  =  'GhostCommunication ' // trim ( FSC % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GC, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GC, Level = 1 )
      end if
    end if
    if ( iT_GPU == 0 ) then
      TimerName  =  'GhostPackUnpack ' // trim ( FSC % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GPU, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT_GPU, Level = 1 )
      end if
    end if
    end associate !-- iT_GC, etc.
      
    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

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

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC

    call Show ( 'Finishing ghost exchange', FSC % IGNORABILITY + 2 )
    call Show ( FSC % Name, 'FieldSet', FSC % IGNORABILITY + 2 )

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


  subroutine UpdateDevice_FS ( FSC, TimerLevelOption )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  FSC % iTimerUpdateDevice )
    if ( iT == 0 ) then
      TimerName  =  'UpdateDevice ' // trim ( FSC % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( FSC % iTimerUpdateDevice )

    call T % Start ( )
    call FSC % Storage % UpdateDevice ( )
    call T % Stop ( )

  end subroutine UpdateDevice_FS


  subroutine UpdateHost_FS ( FSC, TimerLevelOption )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  FSC % iTimerUpdateHost )
    if ( iT == 0 ) then
      TimerName  =  'UpdateHost ' // trim ( FSC % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( FSC % iTimerUpdateHost )

    call T % Start ( )
    call FSC % Storage % UpdateHost ( )
    call T % Stop ( )

  end subroutine UpdateHost_FS


  impure elemental subroutine Finalize ( FSC )

    type ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC

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

    if ( allocated ( FSC % Storage ) ) &
      deallocate ( FSC % Storage )

  end subroutine Finalize


  subroutine AllocateStorage ( FSC )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    if ( .not. allocated ( FSC % Storage ) ) then
      allocate ( FSC % Storage )
      associate ( S  =>  FSC % Storage )
      call S % Initialize &
             ( [ C % nCellsLocal, FSC % nFields ], &
               VariableOption = FSC % Field, &
               VectorOption = FSC % Vector, &
               NameOption = FSC % Name, &
               ClearOption = .true., &
               PinnedOption = FSC % PinnedMemory, &
               UnitOption = FSC % Unit, &
               VectorIndicesOption = FSC % VectorIndices )
      if ( FSC % DeviceMemory ) &
        call S % AllocateDevice ( )
      end associate !-- S
    end if

    end select !-- G

  end subroutine AllocateStorage


  subroutine CloneStorage ( FSC_T, FSC_S )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
      FSC_T
    class ( FieldSet_CGS_Form ), intent ( in ) :: &
      FSC_S

    if ( .not. allocated ( FSC_T % Storage ) ) then
      allocate ( FSC_T % Storage )
      associate ( S  =>  FSC_T % Storage )
      call S % Initialize &
             ( FSC_S % Storage, &
               VectorOption = FSC_T % Vector, &
               NameOption = FSC_T % Name, &
               VectorIndicesOption = FSC_T % VectorIndices, &
               iaSelectedOption = FSC_T % iaSelected )
      end associate !-- S
    end if

  end subroutine CloneStorage


  subroutine StartExchangeFace &
               ( FSC, IncomingFace, OutgoingFace, PH, TagReceive, TagSend )

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
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
    class is ( Chart_GS_Form )

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
               PH % nChunksFrom  *  FSC % nFields )
      call OutgoingFace % Initialize &
             ( Communicator, TagSend ( : nD ), PH % Target, &
               PH % nChunksTo  *  FSC % nFields )
    
      if ( FSC % DevicesCommunicate ) then
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

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
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
    class is ( Chart_GS_Form )

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

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
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
    class is ( Chart_GS_Form )

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
      
      !-- Post Receives

      call IncomingEdge % Initialize &
             ( Communicator, pack ( TagReceive, DimensionMask ), PH % Source, &
               PH % nChunksFrom  *  FSC % nFields )           
      call OutgoingEdge % Initialize &
             ( Communicator, pack ( TagSend, DimensionMask ), PH % Target, &
               PH % nChunksTo  *  FSC % nFields )
      
      if ( FSC % DevicesCommunicate ) then
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

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
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
    class is ( Chart_GS_Form )

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

    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
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
    class is ( Chart_GS_Form )

    associate ( S  =>  FSC % Storage )

    oBuffer = 0
    do iS = 1, S % nVariables
      iF = S % iaSelected ( iS )
      call C % SetFieldPointer ( S % Value ( :, iF ), F )
      call Copy ( F, nSend, oSend, oBuffer, OutgoingMessage % Value, &
                  UseDeviceOption = FSC % DevicesCommunicate )
      oBuffer = oBuffer + product ( nSend )
    end do !-- iS

    end associate !-- S
    end select !-- C
    nullify ( F )

    call T % Stop ( )
    nullify ( T )

  end subroutine LoadMessage


  subroutine StoreMessage ( FSC, IncomingMessage, nReceive, oReceive )
               
    class ( FieldSet_CGS_Form ), intent ( inout ) :: &
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
    class is ( Chart_GS_Form )

    associate ( S  =>  FSC % Storage )

    oBuffer = 0
    do iS = 1, S % nVariables          
      iF = S % iaSelected ( iS )
      call C % SetFieldPointer ( S % Value ( :, iF ), F )
      call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, F, &
                  UseDeviceOption = FSC % DevicesCommunicate )
      oBuffer = oBuffer + product ( nReceive )
    end do !-- iS
    
    end associate !-- S
    end select !-- C
    nullify ( F )

    call T % Stop ( )    
    nullify ( T )

  end subroutine StoreMessage


end module FieldSet_CGS__Form
