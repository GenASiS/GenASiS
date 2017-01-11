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
    type ( VariableGroupForm ), dimension ( : ), allocatable :: &
      ExchangeVariableGroup
    type ( VariableGroup_1D_Form ), allocatable :: &
      VariableGroup_1D
    type ( PortalHeaderForm ), allocatable :: &
      PortalPreviousNext, &
      PortalNextPrevious
    type ( MessageIncoming_1D_R_Form ), allocatable :: &
      Incoming
    type ( MessageOutgoing_1D_R_Form ), allocatable :: &
      Outgoing
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
      StartExchange, &
      FinishExchange

      private :: &
        LoadMessage, &
        StoreMessage

    integer ( KDI ), dimension ( 3 ), private, parameter :: &
      TAG_RECEIVE_PREVIOUS = [ 99, 98, 97 ], &
      TAG_RECEIVE_NEXT     = [ 96, 95, 94 ], &
      TAG_SEND_PREVIOUS    = TAG_RECEIVE_NEXT, &
      TAG_SEND_NEXT        = TAG_RECEIVE_PREVIOUS

contains


  subroutine InitializeBasic &
               ( C, Atlas, IsPeriodic, iChart, SpacingOption, &
                 CoordinateSystemOption, IsDistributedOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, nCellsOption, &
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
      SpacingOption
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
           ( Atlas, IsPeriodic, iChart, SpacingOption, &
             CoordinateSystemOption, IsDistributedOption, &
             CoordinateUnitOption, MinCoordinateOption, &
             MaxCoordinateOption, RatioOption, ScaleOption, nCellsOption, &
             nGhostLayersOption, nDimensionsOption, nEqualOption )

    call SetGhostExchangePortals ( C )

  end subroutine InitializeBasic


  subroutine ExchangeGhostSingle ( C, VG )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    class ( VariableGroupForm ), intent ( inout ) :: &
      VG

    allocate ( C % ExchangeVariableGroup ( 1 ) )
    call C % ExchangeVariableGroup ( 1 ) % Initialize ( VG )

    call ExchangeGhostMultiple ( C, C % ExchangeVariableGroup )
    
  end subroutine ExchangeGhostSingle


  subroutine ExchangeGhostMultiple ( C, VG )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      VG

    call StartExchange &
           ( C, C % PortalPreviousNext, VG, TAG_RECEIVE_PREVIOUS, &
             TAG_SEND_NEXT )
    call FinishExchange ( C, TAG_RECEIVE_PREVIOUS )

    call StartExchange &
           ( C, C % PortalNextPrevious, VG, TAG_RECEIVE_NEXT, &
             TAG_SEND_PREVIOUS )
    call FinishExchange ( C, TAG_RECEIVE_NEXT )

    if ( allocated ( C % ExchangeVariableGroup ) ) &
      deallocate ( C % ExchangeVariableGroup )

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

    if ( allocated ( C % Outgoing ) ) &
      deallocate ( C % Outgoing )
    if ( allocated ( C % Incoming ) ) &
      deallocate ( C % Incoming )
    if ( allocated ( C % PortalNextPrevious ) ) &
      deallocate ( C % PortalNextPrevious )
    if ( allocated ( C % PortalPreviousNext ) ) &
      deallocate ( C % PortalPreviousNext )
    if ( allocated ( C % VariableGroup_1D ) ) &
      deallocate ( C % VariableGroup_1D )
    if ( allocated ( C % ExchangeVariableGroup ) ) &
      deallocate ( C % ExchangeVariableGroup )

    call C % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetGeometryWidthCenter ( C, Width, Center, iD )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( Real_1D_Form ), intent ( inout ) :: &
      Width, &
      Center
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension

    integer ( KDI ) :: &
      iC, &  !-- iCell
      oC  !-- oCell (offset)
    type ( Real_1D_Form ) :: &
      GlobalWidth, &
      GlobalCenter

    !-- Global cell widths and centers

    associate &
      ( iaF => 1 - C % nGhostLayers ( iD ), &
        iaL => C % nCells ( iD ) + C % nGhostLayers ( iD ) )
    call GlobalWidth  % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    call GlobalCenter % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    end associate !-- iaF, etc.

    associate &
      ( Width_1D_G  => GlobalWidth  % Value, &
        Center_1D_G => GlobalCenter % Value, &
        nC => C % nCells ( iD ), &
        nGL => C % nGhostLayers ( iD ) )

    associate ( iaF => 1 - C % nGhostLayers ( iD ) )
    call C % SetGeometryCell ( Width_1D_G, Center_1D_G, nC, nGL, iD, iaF )
    end associate !-- iaF

    !-- Local cell widths and centers

    associate &
      ( iaF => 1 - C % nGhostLayers ( iD ), &
        iaL => C % nCellsBrick ( iD ) + C % nGhostLayers ( iD ) )
    call Width  % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    call Center % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )

    associate &          
      ( Width_1D_L  => Width  % Value, &
        Center_1D_L => Center % Value )

    oC = ( C % iaBrick ( iD ) - 1 ) * C % nCellsBrick ( iD )
    do iC = iaF, iaL
      Width_1D_L  ( iC ) = Width_1D_G  ( oC + iC )
      Center_1D_L ( iC ) = Center_1D_G ( oC + iC )
    end do

    end associate !-- Width_1D_L, etc.
    end associate !-- iaF, etc.

    end associate !-- Width_1D_G, etc.

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
      Source_PN, &
      Source_NP, &
      Target_PN, &
      Target_NP
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
    allocate ( Source_PN ( nD ) )
    allocate ( Source_NP ( nD ) )
    allocate ( Target_PN ( nD ) )
    allocate ( Target_NP ( nD ) )
    
    !-- Face sibling bricks

    do iD = 1, nD
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1

      nCells ( iD ) = C % nGhostLayers ( iD ) &
                        * C % nCellsBrick ( jD ) * C % nCellsBrick ( kD )

      iaB = C % iaBrick

      !-- Previous brick

      iaB ( iD ) = mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      Source_PN ( iD ) = Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_NP ( iD ) = Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- Next brick

      iaB ( iD ) = mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
      Source_NP ( iD ) = Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_PN ( iD ) = Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

    end do !-- iD

    !-- Set portals
    
    allocate ( C % PortalPreviousNext )
    associate ( PPN => C % PortalPreviousNext )
    call PPN % Initialize ( Source_PN, Target_PN, nCells, nCells )
    call PPN % Show ( 'DistributedMesh PortalPreviousNext', &
                      C % IGNORABILITY + 1 )
    end associate !-- PPN

    allocate ( C % PortalNextPrevious )
    associate ( PNP => C % PortalNextPrevious )
    call PNP % Initialize ( Source_NP, Target_NP, nCells, nCells )
    call PNP % Show ( 'DistributedMesh PortalNextPrevious', &
                      C % IGNORABILITY + 1 )
    end associate !-- PNP

    end associate !-- nB
    
  end subroutine SetGhostExchangePortals


  subroutine StartExchange ( C, PH, VG, TagReceive, TagSend )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( PortalHeaderForm ), intent ( in ) :: &
      PH
    type ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      VG
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      TagReceive, &
      TagSend

    integer ( KDI ) :: &
      iD !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend

    allocate ( C % VariableGroup_1D )
    allocate ( C % Incoming )
    allocate ( C % Outgoing )

    associate &
      ( VG_1D  => C % VariableGroup_1D, &
        Communicator => C % Atlas % Communicator, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers, &
        nD  => C % nDimensions )

    call VG_1D % Initialize ( VG )

    !-- Post Receives

    call C % Incoming % Initialize &
           ( Communicator, TagReceive ( : nD ), PH % Source, &
             PH % nChunksFrom * VG_1D % nVariablesTotal )
    call C % Incoming % Receive ( )

    !-- Post Sends

    call C % Outgoing % Initialize &
           ( Communicator, TagSend ( : nD ), PH % Target, &
             PH % nChunksTo * VG_1D % nVariablesTotal )

    do iD = 1, nD

      nSend        = nCB
      nSend ( iD ) = nGL ( iD )

      !-- In setting oSend, note Copy command does not inherit lbound
      if ( all ( TagSend == TAG_SEND_NEXT ) ) then
        oSend        = nGL
        oSend ( iD ) = oSend ( iD ) + nCB ( iD ) - nGL ( iD )
      else if ( all ( TagSend == TAG_SEND_PREVIOUS ) ) then
        oSend = nGL
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'Chart_SLD__Form', 'module', CONSOLE % ERROR )
        call Show ( 'StartExchange', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagSend

      call LoadMessage ( C, VG_1D, nSend, oSend, iD )
      call C % Outgoing % Send ( iD )

    end do !-- iD

    end associate !-- VG_1D, etc.

  end subroutine StartExchange


  subroutine FinishExchange ( C, TagReceive )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
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

    associate &
      ( VG_1D  => C % VariableGroup_1D, &
        nCB => C % nCellsBrick, &
        nGL => C % nGhostLayers )

    !-- Wait for Receives

    do 

      call C % Incoming % Wait ( AllFinished, iD )
      if ( AllFinished ) exit

      nReceive        = nCB
      nReceive ( iD ) = nGL ( iD )

      !-- In setting oReceive, note Copy command does not inherit lbound
      if ( all ( TagReceive == TAG_RECEIVE_PREVIOUS ) ) then
        if ( C % iaBrick ( iD ) == 1 &
             .and. .not. C % IsPeriodic ( iD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) - nGL ( iD )
      else if ( all ( TagReceive == TAG_RECEIVE_NEXT ) ) then
        if ( C % iaBrick ( iD ) == C % nBricks ( iD ) &
             .and. .not. C % IsPeriodic ( iD ) ) &
          cycle
        oReceive        = nGL
        oReceive ( iD ) = oReceive ( iD ) + nCB ( iD )
      else
        call Show ( 'Tags not recognized', CONSOLE % ERROR )
        call Show ( 'Chart_SLD__Form', 'module', CONSOLE % ERROR )
        call Show ( 'FinishExchange', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if !-- TagReceive

      call StoreMessage ( C, VG_1D, nReceive, oReceive, iD )

    end do

    !-- Wait for Sends

    call C % Outgoing % Wait ( )

    !-- Cleanup

    end associate !-- VG_1D, etc.

    deallocate ( C % Outgoing )
    deallocate ( C % Incoming )
    deallocate ( C % VariableGroup_1D )

  end subroutine FinishExchange


  subroutine LoadMessage ( C, VG_1D, nSend, oSend, iM )

    class ( Chart_SLD_Form ), intent ( inout ) :: &
      C
    type ( VariableGroup_1D_Form ), intent ( in ) :: &
      VG_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nSend, &
      oSend
    integer ( KDI ), intent ( in ) :: &
      iM

    integer ( KDI ) :: &
      iG, &  !-- iGroup
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable

    oBuffer = 0
    do iG = 1, VG_1D % nGroups
      do iS = 1, VG_1D % nVariables ( iG )          
        iV = VG_1D % VariableGroup ( iG ) % iaSelected ( iS )
        call C % SetVariablePointer &
               ( VG_1D % VariableGroup ( iG ) % Value ( :, iV ), V ) 
        call Copy ( V, nSend, oSend, oBuffer, &
                    C % Outgoing % Message ( iM ) % Value )
        oBuffer = oBuffer + product ( nSend )
      end do !-- iS
    end do !-- iG

    nullify ( V )

  end subroutine LoadMessage


  subroutine StoreMessage ( C, VG_1D, nReceive, oReceive, iM )

    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    type ( VariableGroup_1D_Form ), intent ( inout ) :: &
      VG_1D
    integer ( KDI ), dimension ( 3 ), intent ( in )  :: &
      nReceive, &
      oReceive
    integer ( KDI ), intent ( in ) :: &
      iM

    integer ( KDI ) :: &
      iG, &  !-- iGroup
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable

    oBuffer = 0
    do iG = 1, VG_1D % nGroups
      do iS = 1, VG_1D % nVariables ( iG )          
        iV = VG_1D % VariableGroup ( iG ) % iaSelected ( iS )
        call C % SetVariablePointer &
               ( VG_1D % VariableGroup ( iG ) % Value ( :, iV ), V ) 
        call Copy ( C % Incoming % Message ( iM ) % Value, &
                    nReceive, oReceive, oBuffer, V )
        oBuffer = oBuffer + product ( nReceive )
      end do !-- iS
    end do !-- iG

    nullify ( V )

  end subroutine StoreMessage


end module Chart_SLD__Form
