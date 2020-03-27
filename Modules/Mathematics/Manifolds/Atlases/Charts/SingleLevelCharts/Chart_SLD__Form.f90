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
      PortalFace_R_L
    type ( MessageIncoming_1D_R_Form ), allocatable :: &
      IncomingFace_L_R, &
      IncomingFace_R_L
    type ( MessageOutgoing_1D_R_Form ), allocatable :: &
      OutgoingFace_L_R, &
      OutgoingFace_R_L
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
      FinishExchangeFace

      private :: &
        LoadMessage, &
        StoreMessage

    integer ( KDI ), dimension ( 3 ), private, parameter :: &
      TAG_RECEIVE_FACE_L = [ 99, 98, 97 ], &
      TAG_RECEIVE_FACE_R = [ 96, 95, 94 ], &
      TAG_SEND_FACE_L    = TAG_RECEIVE_FACE_R, &
      TAG_SEND_FACE_R    = TAG_RECEIVE_FACE_L

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

    !-- Finish faces

    call FinishExchangeFace &
           ( C, C % IncomingFace_L_R, C % OutgoingFace_L_R, TAG_RECEIVE_FACE_L )
    call FinishExchangeFace &
           ( C, C % IncomingFace_R_L, C % OutgoingFace_R_L, TAG_RECEIVE_FACE_R )

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
      iD, jD, kD!, &  !-- iDimension, etc.
    integer ( KDI ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      iaB
    integer ( KDI ), dimension ( : ), allocatable :: &
      nCells, &
      Source_L_R, &
      Source_R_L, &
      Target_L_R, &
      Target_R_L
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

    allocate ( nCells ( nD ) )
    allocate ( Source_L_R ( nD ) )
    allocate ( Source_R_L ( nD ) )
    allocate ( Target_L_R ( nD ) )
    allocate ( Target_R_L ( nD ) )
    
    !-- Face sibling bricks

    do iD = 1, nD

      jD  =  mod ( iD, 3 ) + 1
      kD  =  mod ( jD, 3 ) + 1

      nCells ( iD )  =  C % nGhostLayers ( iD ) &
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

    !-- Set portals
    
    allocate ( C % PortalFace_L_R )
    associate ( PPN => C % PortalFace_L_R )
    call PPN % Initialize ( Source_L_R, Target_L_R, nCells, nCells )
    call PPN % Show ( 'DistributedMesh PortalFace_L_R', &
                      C % IGNORABILITY + 2 )
    end associate !-- PPN

    allocate ( C % PortalFace_R_L )
    associate ( PNP => C % PortalFace_R_L )
    call PNP % Initialize ( Source_R_L, Target_R_L, nCells, nCells )
    call PNP % Show ( 'DistributedMesh PortalFace_R_L', &
                      C % IGNORABILITY + 2 )
    end associate !-- PNP

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

      nSend        = nCB
      nSend ( iD ) = nGL ( iD )

      !-- In setting oSend, note Copy command does not inherit lbound
      if ( all ( TagSend == TAG_SEND_FACE_R ) ) then
        oSend        = nGL
        oSend ( iD ) = oSend ( iD ) + nCB ( iD ) - nGL ( iD )
      else if ( all ( TagSend == TAG_SEND_FACE_L ) ) then
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
      if ( all ( TagReceive == TAG_RECEIVE_FACE_L ) ) then
        if ( C % iaBrick ( iD ) == 1 &
             .and. .not. C % IsPeriodic ( iD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) - nGL ( iD )
      else if ( all ( TagReceive == TAG_RECEIVE_FACE_R ) ) then
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
