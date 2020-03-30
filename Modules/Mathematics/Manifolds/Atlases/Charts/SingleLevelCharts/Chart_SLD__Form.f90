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
    type ( StorageForm ), dimension ( : ), allocatable :: &
      ExchangeStorage
    type ( Storage_1D_Form ), allocatable :: &
      Storage_1D
    type ( PortalHeaderForm ), allocatable :: &
      PortalFace_L_R, &
      PortalFace_R_L, &
      PortalEdge_LL_RR, &
      PortalEdge_RR_LL, &
      PortalEdge_LR_RL, &
      PortalEdge_RL_LR
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
      InitializeBasic
    procedure, public, pass :: &
      ExchangeGhostSingle
    procedure, public, pass :: &
      ExchangeGhostMultiple
    generic :: &
      ExchangeGhostData => ExchangeGhostSingle, ExchangeGhostMultiple
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
                 ScaleOption, nCellsOption, nGhostLayersOption, &
                 nDimensionsOption, nEqualOption )

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
             RatioOption, ScaleOption, nCellsOption, nGhostLayersOption, &
             nDimensionsOption, nEqualOption )

    call SetGhostExchangePortals ( C )

  end subroutine InitializeBasic


  subroutine ExchangeGhostSingle ( C, S )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( StorageForm ), intent ( inout ) :: &
      S

    allocate ( C % ExchangeStorage ( 1 ) )
    call C % ExchangeStorage ( 1 ) % Initialize ( S )

    call ExchangeGhostMultiple ( C, C % ExchangeStorage )
    
    deallocate ( C % ExchangeStorage )

  end subroutine ExchangeGhostSingle


  subroutine ExchangeGhostMultiple ( C, S )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( StorageForm ), dimension ( : ), intent ( in ) :: &
      S

    allocate ( C % Storage_1D )
    call C % Storage_1D % Initialize ( S )

    !-- Start faces
    call StartExchangeFace &
           ( C, C % PortalFace_L_R, TAG_RECEIVE_FACE_L, TAG_SEND_FACE_R, &
             C % IncomingFace_L_R, C % OutgoingFace_L_R )
    call StartExchangeFace &
           ( C, C % PortalFace_R_L, TAG_RECEIVE_FACE_R, TAG_SEND_FACE_L, &
             C % IncomingFace_R_L, C % OutgoingFace_R_L )

    !-- Start edges
    call StartExchangeEdge &
           ( C, C % PortalEdge_LL_RR, TAG_RECEIVE_EDGE_LL, TAG_SEND_EDGE_RR, &
             C % IncomingEdge_LL_RR, C % OutgoingEdge_LL_RR )
    call StartExchangeEdge &
           ( C, C % PortalEdge_RR_LL, TAG_RECEIVE_EDGE_RR, TAG_SEND_EDGE_LL, &
             C % IncomingEdge_RR_LL, C % OutgoingEdge_RR_LL )
    call StartExchangeEdge &
           ( C, C % PortalEdge_LR_RL, TAG_RECEIVE_EDGE_LR, TAG_SEND_EDGE_RL, &
             C % IncomingEdge_LR_RL, C % OutgoingEdge_LR_RL )
    call StartExchangeEdge &
           ( C, C % PortalEdge_RL_LR, TAG_RECEIVE_EDGE_RL, TAG_SEND_EDGE_LR, &
             C % IncomingEdge_RL_LR, C % OutgoingEdge_RL_LR )

    !-- Finish faces
    call FinishExchangeFace &
           ( C, C % IncomingFace_L_R, C % OutgoingFace_L_R, TAG_RECEIVE_FACE_L )
    call FinishExchangeFace &
           ( C, C % IncomingFace_R_L, C % OutgoingFace_R_L, TAG_RECEIVE_FACE_R )

    !-- Finish edges
    call FinishExchangeEdge &
           ( C, C % IncomingEdge_LL_RR, C % OutgoingEdge_LL_RR, &
             TAG_RECEIVE_EDGE_LL )
    call FinishExchangeEdge &
           ( C, C % IncomingEdge_RR_LL, C % OutgoingEdge_RR_LL, &
             TAG_RECEIVE_EDGE_RR )
    call FinishExchangeEdge &
           ( C, C % IncomingEdge_LR_RL, C % OutgoingEdge_LR_RL, &
             TAG_RECEIVE_EDGE_LR )
    call FinishExchangeEdge &
           ( C, C % IncomingEdge_RL_LR, C % OutgoingEdge_RL_LR, &
             TAG_RECEIVE_EDGE_RL )

    deallocate ( C % Storage_1D )

  end subroutine ExchangeGhostMultiple


  subroutine CopyBoundary ( Value, C, iDimension, iConnection )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
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
           ( Value, C % nCellsBrick, iDimension, iConnection )

  end subroutine CopyBoundary


  subroutine ReverseBoundary ( Value, C, iDimension, iConnection )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
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
           ( Value, C % nCellsBrick, iDimension, iConnection )

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
               ( C, PH, TagReceive, TagSend, IncomingFace, OutgoingFace )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend
    type ( MessageIncoming_1D_R_Form ), intent ( out ), allocatable :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( out ), allocatable :: &
      OutgoingFace

    integer ( KDI ) :: &
      iD !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend

    allocate ( IncomingFace )
    allocate ( OutgoingFace )

    associate &
      ( S_1D  => C % Storage_1D, &
        Communicator => C % Atlas % Communicator, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        nD  => C % nDimensions )

    !-- Post Receives

    call IncomingFace % Initialize &
           ( Communicator, TagReceive ( : nD ), PH % Source, &
             PH % nChunksFrom * S_1D % nVariablesTotal )
    call IncomingFace % Receive ( )

    !-- Post Sends

    call OutgoingFace % Initialize &
           ( Communicator, TagSend ( : nD ), PH % Target, &
             PH % nChunksTo * S_1D % nVariablesTotal )

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

      call LoadMessage ( C, OutgoingFace % Message ( iD ), S_1D, nSend, oSend )
      call OutgoingFace % Send ( iD )

    end do !-- iD

    end associate !-- S_1D, etc.

  end subroutine StartExchangeFace


  subroutine FinishExchangeFace ( C, IncomingFace, OutgoingFace, TagReceive )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingFace
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
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

    associate &
      ( S_1D  => C % Storage_1D, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers )

    !-- Wait for Receives

    do 

      call IncomingFace % Wait ( AllFinished, iD )
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
             ( C, S_1D, IncomingFace % Message ( iD ), nReceive, oReceive )

    end do

    !-- Wait for Sends

    call OutgoingFace % Wait ( )

    !-- Cleanup

    end associate !-- S_1D, etc.

    deallocate ( OutgoingFace )
    deallocate ( IncomingFace )

  end subroutine FinishExchangeFace


  subroutine StartExchangeEdge &
               ( C, PH, TagReceive, TagSend, IncomingEdge, OutgoingEdge )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend
    type ( MessageIncoming_1D_R_Form ), intent ( out ), allocatable :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( out ), allocatable :: &
      OutgoingEdge

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension, etc.
      kM  !-- kMessage
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend
    logical ( KDL ), dimension ( 3 ) :: &
      DimensionMask

    allocate ( IncomingEdge )
    allocate ( OutgoingEdge )

    associate &
      ( S_1D  => C % Storage_1D, &
        Communicator => C % Atlas % Communicator, &
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

    !-- Post Receives

    call IncomingEdge % Initialize &
           ( Communicator, pack ( TagReceive, DimensionMask ), PH % Source, &
             PH % nChunksFrom * S_1D % nVariablesTotal )
    call IncomingEdge % Receive ( )

    !-- Post Sends

    call OutgoingEdge % Initialize &
           ( Communicator, pack ( TagSend, DimensionMask ), PH % Target, &
             PH % nChunksTo * S_1D % nVariablesTotal )

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

      call LoadMessage ( C, OutgoingEdge % Message ( kM ), S_1D, nSend, oSend )
      call OutgoingEdge % Send ( kM )

    end do !-- kD

    end associate !-- S_1D, etc.

  end subroutine StartExchangeEdge


  subroutine FinishExchangeEdge ( C, IncomingEdge, OutgoingEdge, TagReceive )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( MessageIncoming_1D_R_Form ), intent ( inout ), allocatable :: &
      IncomingEdge
    type ( MessageOutgoing_1D_R_Form ), intent ( inout ), allocatable :: &
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

    associate &
      ( S_1D => C % Storage_1D, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        nD  => C % nDimensions )

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
             ( C, S_1D, IncomingEdge % Message ( kM ), nReceive, oReceive )

    end do

    !-- Wait for Sends

    call OutgoingEdge % Wait ( )

    !-- Cleanup

    end associate !-- S_1D, etc.

    deallocate ( OutgoingEdge )
    deallocate ( IncomingEdge )

  end subroutine FinishExchangeEdge


  subroutine LoadMessage ( C, OutgoingMessage, S_1D, nSend, oSend )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( MessageOutgoing_R_Form ), intent ( in ) :: &
      OutgoingMessage
    type ( Storage_1D_Form ), intent ( in ) :: &
      S_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nSend, &
      oSend

    integer ( KDI ) :: &
      iStrg, &  !-- iStorage
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable

    oBuffer = 0
    do iStrg = 1, S_1D % nStorages
      do iS = 1, S_1D % nVariables ( iStrg )          
        iV = S_1D % Storage ( iStrg ) % iaSelected ( iS )
        call C % SetVariablePointer &
               ( S_1D % Storage ( iStrg ) % Value ( :, iV ), V ) 
        call Copy ( V, nSend, oSend, oBuffer, OutgoingMessage % Value )
        oBuffer = oBuffer + product ( nSend )
      end do !-- iS
    end do !-- iStrg

    nullify ( V )

  end subroutine LoadMessage


  subroutine StoreMessage ( C, S_1D, IncomingMessage, nReceive, oReceive )

    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    type ( Storage_1D_Form ), intent ( inout ) :: &
      S_1D
    type ( MessageIncoming_R_Form ), intent ( in ) :: &
      IncomingMessage
    integer ( KDI ), dimension ( 3 ), intent ( in )  :: &
      nReceive, &
      oReceive

    integer ( KDI ) :: &
      iStrg, &  !-- iStorage
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable

    oBuffer = 0
    do iStrg = 1, S_1D % nStorages
      do iS = 1, S_1D % nVariables ( iStrg )          
        iV = S_1D % Storage ( iStrg ) % iaSelected ( iS )
        call C % SetVariablePointer &
               ( S_1D % Storage ( iStrg ) % Value ( :, iV ), V ) 
        call Copy ( IncomingMessage % Value, nReceive, oReceive, oBuffer, V )
        oBuffer = oBuffer + product ( nReceive )
      end do !-- iS
    end do !-- iStrg

    nullify ( V )

  end subroutine StoreMessage


end module Chart_SLD__Form
