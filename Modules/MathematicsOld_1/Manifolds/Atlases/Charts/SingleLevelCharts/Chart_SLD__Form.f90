!-- Chart_SLD represents a distributed single-level chart.

module Chart_SLD__Form

  !-- Chart_SingleLevelDistributed_Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SL__Template

  implicit none
  private

  type, public, extends ( Chart_SL_Template ) :: Chart_SLD_Form
    logical ( KDL ) :: &
      ExchangeGhostUseDevice
    type ( PortalHeaderForm ), allocatable :: &
      PortalFace_L_R, &
      PortalFace_R_L, &
      PortalEdge_LL_RR, &
      PortalEdge_RR_LL, &
      PortalEdge_LR_RL, &
      PortalEdge_RL_LR
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, public, pass :: &
      ExchangeGhostData
    procedure, public, pass ( C ) :: &
      CopyBoundary
    procedure, public, pass ( C ) :: &
      ReverseBoundary
    final :: &
      Finalize
    procedure, public, pass :: &
      SetGeometryWidthCenter
  end type Chart_SLD_Form

    private :: &
      SetGhostExchangePortals, &
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


  subroutine InitializeBasic &
               ( C, Atlas, IsPeriodic, iChart, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 IsDistributedOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nBricksOption, &
                 nGhostLayersOption, nDimensionsOption, nEqualOption )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsPeriodic
    integer ( KDI ), intent ( in ) :: &
      iChart
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    logical ( KDL ), intent ( in ), optional :: &
      IsDistributedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nBricksOption, &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    if ( .not. associated ( C % Type ) ) then
      allocate ( C % Type )
      C % Type = 'a Chart_SLD'
    end if

    call C % ChartHeader_SL_Form % Initialize &
           ( Atlas, IsPeriodic, iChart, SpacingOption, CoordinateLabelOption, &
             CoordinateSystemOption, IsDistributedOption, &
             CoordinateUnitOption, MinCoordinateOption, MaxCoordinateOption, &
             RatioOption, ScaleOption, nCellsOption, nBricksOption, &
             nGhostLayersOption, nDimensionsOption, nEqualOption )
    
    C % ExchangeGhostUseDevice = .true.
    
    call PROGRAM_HEADER % GetParameter &
           ( C % ExchangeGhostUseDevice, 'ExchangeGhostUseDevice', &
             IgnorabilityOption = CONSOLE % INFO_2 )

    call SetGhostExchangePortals ( C )

  end subroutine InitializeBasic


  subroutine ExchangeGhostData ( C, F, UseDeviceOption )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( Field_CSL_Template ), intent ( inout ) :: &
      F
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    associate ( S => F % Field )
    
    !-- Start faces
    call StartExchangeFace &
           ( C, S, F % IncomingFace_L_R, F % OutgoingFace_L_R, &
             C % PortalFace_L_R, UseDevice, &
             TAG_RECEIVE_FACE_L, TAG_SEND_FACE_R )
    call StartExchangeFace &
           ( C, S, F % IncomingFace_R_L, F % OutgoingFace_R_L, &
             C % PortalFace_R_L, UseDevice, &
             TAG_RECEIVE_FACE_R, TAG_SEND_FACE_L )

    !-- Start edges
    call StartExchangeEdge &
           ( C, S, F % IncomingEdge_LL_RR, F % OutgoingEdge_LL_RR, &
             C % PortalEdge_LL_RR, UseDevice, &
             TAG_RECEIVE_EDGE_LL, TAG_SEND_EDGE_RR )
    call StartExchangeEdge &
           ( C, S, F % IncomingEdge_RR_LL, F % OutgoingEdge_RR_LL, &
             C % PortalEdge_RR_LL, UseDevice, &
             TAG_RECEIVE_EDGE_RR, TAG_SEND_EDGE_LL )
    call StartExchangeEdge &
           ( C, S, F % IncomingEdge_LR_RL, F % OutgoingEdge_LR_RL, &
             C % PortalEdge_LR_RL, UseDevice, &
             TAG_RECEIVE_EDGE_LR, TAG_SEND_EDGE_RL )
    call StartExchangeEdge &
           ( C, S, F % IncomingEdge_RL_LR, F % OutgoingEdge_RL_LR, &
             C % PortalEdge_RL_LR, UseDevice, &
             TAG_RECEIVE_EDGE_RL, TAG_SEND_EDGE_LR )

    !-- Finish faces
    call FinishExchangeFace &
           ( C, S, F % IncomingFace_L_R, F % OutgoingFace_L_R, &
             UseDevice, TAG_RECEIVE_FACE_L )
    call FinishExchangeFace &
           ( C, S, F % IncomingFace_R_L, F % OutgoingFace_R_L, &
             UseDevice, TAG_RECEIVE_FACE_R )

    !-- Finish edges
    call FinishExchangeEdge &
           ( C, S, F % IncomingEdge_LL_RR, F % OutgoingEdge_LL_RR, &
             UseDevice, TAG_RECEIVE_EDGE_LL )
    call FinishExchangeEdge &
           ( C, S, F % IncomingEdge_RR_LL, F % OutgoingEdge_RR_LL, &
             UseDevice, TAG_RECEIVE_EDGE_RR )
    call FinishExchangeEdge &
           ( C, S, F % IncomingEdge_LR_RL, F % OutgoingEdge_LR_RL, &
             UseDevice, TAG_RECEIVE_EDGE_LR )
    call FinishExchangeEdge &
           ( C, S, F % IncomingEdge_RL_LR, F % OutgoingEdge_RL_LR, &
             UseDevice, TAG_RECEIVE_EDGE_RL )
    
    end associate

  end subroutine ExchangeGhostData


  subroutine CopyBoundary ( F, C, iField, iDimension, iConnection )

    class ( StorageForm ), intent ( inout ) :: &
      F
    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iField, &
      iDimension, &
      iConnection

    associate &
      ( Connectivity => C % Atlas % Connectivity, &
        iD => iDimension, &
        iC => iConnection )
    if ( iC == Connectivity % iaInner ( iD ) &
         .and. C % iaBrick ( iD ) /= 1 ) return
    if ( iC == Connectivity % iaOuter ( iD ) &
         .and. C % iaBrick ( iD ) /= C % nBricks ( iD ) ) return
    end associate !-- Connectivity, etc.

    call C % CopyBoundaryTemplate &
           ( F, C % nCellsBrick, iField, iDimension, iConnection )

  end subroutine CopyBoundary


  subroutine ReverseBoundary ( F, C, iField, iDimension, iConnection )

    class ( StorageForm ), intent ( inout ) :: &
      F
    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iField, &
      iDimension, &
      iConnection

    associate &
      ( Connectivity => C % Atlas % Connectivity, &
        iD => iDimension, &
        iC => iConnection )
    if ( iC == Connectivity % iaInner ( iD ) &
         .and. C % iaBrick ( iD ) /= 1 ) return
    if ( iC == Connectivity % iaOuter ( iD ) &
         .and. C % iaBrick ( iD ) /= C % nBricks ( iD ) ) return
    end associate !-- Connectivity, etc.

    call C % ReverseBoundaryTemplate &
           ( F, C % nCellsBrick, iField, iDimension, iConnection )

  end subroutine ReverseBoundary


  impure elemental subroutine Finalize ( C )

    type ( Chart_SLD_Form ), intent ( inout ) :: &
      C

    if ( allocated ( C % PortalFace_R_L ) ) &
      deallocate ( C % PortalFace_R_L )
    if ( allocated ( C % PortalFace_L_R ) ) &
      deallocate ( C % PortalFace_L_R )

    call C % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetGeometryWidthCenter &
               ( C, Center, Width_L, Width_R, iD, EdgeValueOption )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( Real_1D_Form ), intent ( inout ) :: &
      Center, &
      Width_L, &
      Width_R
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      EdgeValueOption

    integer ( KDI ) :: &
      iC, &  !-- iCell
      oC  !-- oCell (offset)

    !-- Global cell widths and centers

    associate &
      ( nC => C % nCells ( iD ), &
        nGL => C % nGhostLayers ( iD ) )
    call C % SetGeometryCell ( nC, nGL, iD, EdgeValueOption )
    end associate !-- nC, etc.

    !-- Local cell widths and centers

    associate &
      ( iaF => 1 - C % nGhostLayers ( iD ), &
        iaL => C % nCellsBrick ( iD ) + C % nGhostLayers ( iD ) )

    call Center % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    call Width_L % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    call Width_R % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )

    oC = ( C % iaBrick ( iD ) - 1 ) * C % nCellsBrick ( iD )
    do iC = iaF, iaL
      Center % Value ( iC )  = C % Center ( iD ) % Value ( oC + iC )
      Width_L % Value ( iC ) = C % WidthLeft ( iD ) % Value ( oC + iC )
      Width_R % Value ( iC ) = C % WidthRight ( iD ) % Value ( oC + iC )
    end do

    end associate !-- iaF, etc.

  end subroutine SetGeometryWidthCenter


  subroutine SetGhostExchangePortals ( C )

    class ( Chart_SLD_Form ) , intent ( inout )  :: &
      C

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iB, jB, kB, &  !-- iBrick, jBrick, kBrick
      iD, jD, kD     !-- iDimension, etc.
    integer ( KDI ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      iaB
    integer ( KDI ), dimension ( : ), allocatable :: &
      !-- Faces
      nCellsFace, &
      Source_L_R, &
      Source_R_L, &
      Target_L_R, &
      Target_R_L, &
      !-- Edges
      nCellsEdge, &
      Source_LL_RR, &
      Source_RR_LL, &
      Source_LR_RL, &
      Source_RL_LR, &
      Target_LL_RR, &
      Target_RR_LL, &
      Target_LR_RL, &
      Target_RL_LR
    integer ( KDI ), dimension ( :, :, : ), allocatable :: &
      Process

    associate &
      ( nB  => C % nBricks, &
        nD  => C % nDimensions ) 

    allocate ( Process ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

    iP = 0
    do kB = 1, nB ( 3 )
      do jB = 1, nB ( 2 )
        do iB = 1, nB ( 1 )
          Process ( iB, jB, kB ) = iP
          iP = iP + 1
        end do !-- iB
      end do !-- jB
    end do !-- kB


    !-- Face sibling bricks

    allocate ( nCellsFace ( nD ) )
    allocate ( Source_L_R ( nD ) )
    allocate ( Source_R_L ( nD ) )
    allocate ( Target_L_R ( nD ) )
    allocate ( Target_R_L ( nD ) )
    nCellsFace =  0
    Source_L_R = -1
    Source_R_L = -1
    Target_L_R = -1
    Target_R_L = -1
    
    do iD = 1, nD

      jD  =  mod ( iD, 3 ) + 1
      kD  =  mod ( jD, 3 ) + 1

      nCellsFace ( iD )  =  C % nGhostLayers ( iD ) &
                            * C % nCellsBrick ( jD ) * C % nCellsBrick ( kD )

      iaB  =  C % iaBrick

      !-- Left brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      Source_L_R ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_R_L ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- Right brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
      Source_R_L ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_L_R ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

    end do !-- iD


    !-- Edge sibling bricks

    allocate ( nCellsEdge ( 3 ) )
    allocate ( Source_LL_RR ( 3 ) )
    allocate ( Source_RR_LL ( 3 ) )
    allocate ( Source_LR_RL ( 3 ) )
    allocate ( Source_RL_LR ( 3 ) )
    allocate ( Target_LL_RR ( 3 ) )
    allocate ( Target_RR_LL ( 3 ) )
    allocate ( Target_LR_RL ( 3 ) )
    allocate ( Target_RL_LR ( 3 ) )
    nCellsEdge =  0
    Source_LL_RR = -1
    Source_RR_LL = -1
    Source_LR_RL = -1
    Source_RL_LR = -1
    Target_LL_RR = -1
    Target_RR_LL = -1
    Target_LR_RL = -1
    Target_RL_LR = -1

    do kD = 3, 1, -1

      iD  =  mod ( kD, 3 ) + 1
      jD  =  mod ( iD, 3 ) + 1

      if ( iD > nD .or. jD > nD ) &
        cycle

      nCellsEdge ( kD )  &
        =  C % nGhostLayers ( iD )  *  C % nGhostLayers ( jD )  &
           *  C % nCellsBrick ( kD )

      iaB  =  C % iaBrick

      !-- LeftLeft brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ) - 1 + nB ( jD ) - 1, nB ( jD ) ) + 1
      Source_LL_RR ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_RR_LL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- RightRight brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ), nB ( jD ) ) + 1
      Source_RR_LL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_LL_RR ( KD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- LeftRight brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ), nB ( jD ) ) + 1
      Source_LR_RL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_RL_LR ( KD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- RightLeft brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ) - 1 + nB ( jD ) - 1, nB ( jD ) ) + 1
      Source_RL_LR ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_LR_RL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

    end do !-- kD


    !-- Set face portals
    
    allocate ( C % PortalFace_L_R )
    associate ( PFLR => C % PortalFace_L_R )
    call PFLR % Initialize ( Source_L_R, Target_L_R, nCellsFace, nCellsFace )
    call PFLR % Show ( 'Chart_SLD PortalFace_L_R', &
                       C % IGNORABILITY + 2 )
    end associate !-- PFLR

    allocate ( C % PortalFace_R_L )
    associate ( PFRL => C % PortalFace_R_L )
    call PFRL % Initialize ( Source_R_L, Target_R_L, nCellsFace, nCellsFace )
    call PFRL % Show ( 'Chart_SLD PortalFace_R_L', &
                       C % IGNORABILITY + 2 )
    end associate !-- PFRL


    !-- Set edge portals
    
    allocate ( C % PortalEdge_LL_RR )
    associate ( PELLRR => C % PortalEdge_LL_RR )
    call PELLRR % Initialize &
           ( pack ( Source_LL_RR, Source_LL_RR >= 0 ), &
             pack ( Target_LL_RR, Target_LL_RR >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    call PELLRR % Show ( 'Chart_SLD PortalEdge_LL_RR', &
                         C % IGNORABILITY + 2 )
    end associate !-- PELLRR

    allocate ( C % PortalEdge_RR_LL )
    associate ( PERRLL => C % PortalEdge_RR_LL )
    call PERRLL % Initialize &
           ( pack ( Source_RR_LL, Source_RR_LL >= 0 ), &
             pack ( Target_RR_LL, Target_RR_LL >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    call PERRLL % Show ( 'Chart_SLD PortalEdge_RR_LL', &
                         C % IGNORABILITY + 2 )
    end associate !-- PERRLL

    allocate ( C % PortalEdge_LR_RL )
    associate ( PELRRL => C % PortalEdge_LR_RL )
    call PELRRL % Initialize &
           ( pack ( Source_LR_RL, Source_LR_RL >= 0 ), &
             pack ( Target_LR_RL, Target_LR_RL >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    call PELRRL % Show ( 'Chart_SLD PortalEdge_LR_RL', &
                         C % IGNORABILITY + 2 )
    end associate !-- PELRRL

    allocate ( C % PortalEdge_RL_LR )
    associate ( PERLLR => C % PortalEdge_RL_LR )
    call PERLLR % Initialize &
           ( pack ( Source_RL_LR, Source_RL_LR >= 0 ), &
             pack ( Target_RL_LR, Target_RL_LR >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    call PERLLR % Show ( 'Chart_SLD PortalEdge_RL_LR', &
                         C % IGNORABILITY + 2 )
    end associate !-- PELRRL


    !-- Cleanup

    end associate !-- nB, etc.
    
  end subroutine SetGhostExchangePortals


  subroutine StartExchangeFace &
               ( C, S, IncomingFace, OutgoingFace, PH, UseDevice, &
                 TagReceive, TagSend )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( StorageForm ), intent ( inout ) :: &
      S
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingFace
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    logical ( KDL ), intent ( in ) :: &
      UseDevice
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend
    
    integer ( KDI ) :: &
      iD !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend
    type ( TimerForm ), pointer :: &
      T 
      
    T => PROGRAM_HEADER % TimerPointer ( C % iTimerGhostCommunication )
    
    associate &
      ( Communicator => C % Atlas % Communicator, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        nD  => C % nDimensions )

    !-- Allocate on first use
    if ( .not. allocated ( IncomingFace ) &
         .and. .not. allocated ( OutgoingFace ) ) then
    
      allocate ( IncomingFace )
      allocate ( OutgoingFace )

      call IncomingFace % Initialize &
             ( Communicator, TagReceive ( : nD ), PH % Source, &
               PH % nChunksFrom * S % nVariables )
               
      call OutgoingFace % Initialize &
           ( Communicator, TagSend ( : nD ), PH % Target, &
             PH % nChunksTo * S % nVariables )
    
      if ( UseDevice ) then
        call IncomingFace % AllocateDevice ( )
        call OutgoingFace % AllocateDevice ( )
      end if 
    
    end if
    
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
        oSend        = nGL
        oSend ( iD ) = oSend ( iD ) + nCB ( iD ) - nGL ( iD )
      else if ( TagSend ( iD )  ==  TAG_SEND_FACE_L ( iD ) ) then
        oSend = nGL
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'Chart_SLD__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartExchangeFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      call LoadMessage &
             ( C, OutgoingFace % Message ( iD ), S, UseDevice, &
               nSend, oSend )
      
      call T % Start ( )
      call OutgoingFace % Send ( iD )
      call T % Stop ( )

    end do !-- iD

    end associate !-- Communicator, etc.
    
    nullify ( T )

  end subroutine StartExchangeFace


  subroutine FinishExchangeFace &
               ( C, S, IncomingFace, OutgoingFace, UseDevice, TagReceive )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( StorageForm ), intent ( inout ) :: &
      S
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingFace
    logical ( KDL ), intent ( in ) :: &
      UseDevice
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
      
    T => PROGRAM_HEADER % TimerPointer ( C % iTimerGhostCommunication )

    associate &
      ( nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers )

    !-- Wait for Receives

    do 

      call T % Start ( )
      call IncomingFace % Wait ( AllFinished, iD )
      call T % Stop ( )
      
      if ( AllFinished ) exit

      nReceive        = nCB
      nReceive ( iD ) = nGL ( iD )

      !-- In setting oReceive, note Copy command does not inherit lbound
      if ( TagReceive ( iD )  == TAG_RECEIVE_FACE_L ( iD ) ) then
        if ( C % iaBrick ( iD ) == 1 &
             .and. .not. C % IsPeriodic ( iD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) - nGL ( iD )
      else if ( TagReceive ( iD ) == TAG_RECEIVE_FACE_R ( iD ) ) then
        if ( C % iaBrick ( iD ) == C % nBricks ( iD ) &
             .and. .not. C % IsPeriodic ( iD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) + nCB ( iD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'Chart_SLD__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchangeFace', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage &
             ( C, S, IncomingFace % Message ( iD ), UseDevice, &
               nReceive, oReceive )

    end do

    !-- Wait for Sends
    call T % Start ( )
    call OutgoingFace % Wait ( )
    call T % Stop ( )

    !-- Cleanup

    end associate !-- nCB etc.
    
    nullify ( T )

    !deallocate ( OutgoingFace )
    !deallocate ( IncomingFace )

  end subroutine FinishExchangeFace


  subroutine StartExchangeEdge &
               ( C, S, IncomingEdge, OutgoingEdge, PH, UseDevice, &
                 TagReceive, TagSend )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( StorageForm ), intent ( inout ) :: &
      S
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingEdge
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    logical ( KDL ), intent ( in ) :: &
      UseDevice
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
      
    T => PROGRAM_HEADER % TimerPointer ( C % iTimerGhostCommunication )

    associate &
      ( Communicator => C % Atlas % Communicator, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        nD  => C % nDimensions )
        
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
               PH % nChunksFrom * S % nVariables )
               
      call OutgoingEdge % Initialize &
             ( Communicator, pack ( TagSend, DimensionMask ), PH % Target, &
               PH % nChunksTo * S % nVariables )
      
      if ( UseDevice ) then
        call IncomingEdge % AllocateDevice ( )
        call OutgoingEdge % AllocateDevice ( )
      end if
    
    end if
      
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
        call Show ( 'Chart_SLD__Form', 'module', CONSOLE % ERROR )
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
             ( C, OutgoingEdge % Message ( kM ), S, UseDevice, &
               nSend, oSend )
      
      call T % Start ( )
      call OutgoingEdge % Send ( kM )
      call T % Stop ( )

    end do !-- kD

    end associate !-- Communicator, etc.
    
    nullify ( T )

  end subroutine StartExchangeEdge


  subroutine FinishExchangeEdge &
               ( C, S, IncomingEdge, OutgoingEdge, UseDevice, TagReceive )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( StorageForm ), intent ( inout ) :: &
      S
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
      OutgoingEdge
    logical ( KDL ), intent ( in ) :: &
      UseDevice
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
      
    T => PROGRAM_HEADER % TimerPointer ( C % iTimerGhostCommunication )

    associate &
      ( nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        nD  => C % nDimensions )

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
        if ( C % iaBrick ( iD ) == 1 &
             .and. C % iaBrick ( jD ) == 1 &
             .and. .not. C % IsPeriodic ( iD ) &
             .and. .not. C % IsPeriodic ( jD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) - nGL ( iD )
        oReceive ( jD ) = oReceive ( jD ) - nGL ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RR ( kM ) ) then
        if ( C % iaBrick ( iD ) == C % nBricks ( iD ) &
             .and. C % iaBrick ( jD ) == C % nBricks ( jD ) &
             .and. .not. C % IsPeriodic ( iD ) &
             .and. .not. C % IsPeriodic ( jD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) + nCB ( iD )
        oReceive ( jD ) = oReceive ( jD ) + nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_LR ( kM ) ) then
        if ( C % iaBrick ( iD ) == 1 &
             .and. C % iaBrick ( jD ) == C % nBricks ( jD ) &
             .and. .not. C % IsPeriodic ( iD ) &
             .and. .not. C % IsPeriodic ( jD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) - nGL ( iD )
        oReceive ( jD ) = oReceive ( jD ) + nCB ( jD )
      else if ( TagReceive ( kM )  ==  TAG_RECEIVE_EDGE_RL ( kM ) ) then
        if ( C % iaBrick ( iD ) == C % nBricks ( iD ) &
             .and. C % iaBrick ( jD ) == 1 &
             .and. .not. C % IsPeriodic ( iD ) &
             .and. .not. C % IsPeriodic ( jD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) + nCB ( iD )
        oReceive ( jD ) = oReceive ( jD ) - nGL ( jD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'Chart_SLD__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchangeEdge', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage &
             ( C, S, IncomingEdge % Message ( kM ), UseDevice, &
               nReceive, oReceive )

    end do

    !-- Wait for Sends
    call T % Start ( )
    call OutgoingEdge % Wait ( )
    call T % Stop ( )

    !-- Cleanup

    end associate !-- nCB, etc.
    
    nullify ( T )

    !deallocate ( OutgoingEdge )
    !deallocate ( IncomingEdge )

  end subroutine FinishExchangeEdge


  subroutine LoadMessage ( C, OutgoingMessage, S, UseDevice, nSend, oSend )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( MessageOutgoing_R_Form ), intent ( in ) :: &
      OutgoingMessage
    type ( StorageForm ), intent ( in ) :: &
      S
    logical ( KDL ), intent ( in ) :: &
      UseDevice
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nSend, &
      oSend

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable
    type ( TimerForm ), pointer :: &
      T
    T => PROGRAM_HEADER % TimerPointer ( C % iTimerGhostPackUnpack )
    
    call T % Start ( )
    
    oBuffer = 0
    do iS = 1, S % nVariables
      iV = S % iaSelected ( iS )
      call C % SetVariablePointer ( S % Value ( :, iV ), V )
      call Copy ( V, nSend, oSend, oBuffer, OutgoingMessage % Value, &
                  UseDeviceOption = UseDevice )
      oBuffer = oBuffer + product ( nSend )
    end do !-- iS
    
    call T % Stop ( )
    
    nullify ( T )

    nullify ( V )

  end subroutine LoadMessage


  subroutine StoreMessage &
               ( C, S, IncomingMessage, UseDevice, nReceive, oReceive )

    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    type ( StorageForm ), intent ( inout ) :: &
      S
    type ( MessageIncoming_R_Form ), intent ( in ) :: &
      IncomingMessage
    logical ( KDL ), intent ( in ) :: &
      UseDevice
    integer ( KDI ), dimension ( 3 ), intent ( in )  :: &
      nReceive, &
      oReceive

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable
    type ( TimerForm ), pointer :: &
      T
    
    T => PROGRAM_HEADER % TimerPointer ( C % iTimerGhostPackUnpack )
    
    call T % Start ( )
    
    oBuffer = 0
    do iS = 1, S % nVariables          
      iV = S % iaSelected ( iS )
      call C % SetVariablePointer ( S % Value ( :, iV ), V )
      call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, V, &
                  UseDeviceOption = UseDevice )
      oBuffer = oBuffer + product ( nReceive )
    end do !-- iS
    
    call T % Stop ( )
    
    nullify ( T )

    nullify ( V )

  end subroutine StoreMessage


end module Chart_SLD__Form
