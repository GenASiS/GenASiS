module DistributedMesh_Form

  use Basics

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    MAX_N_DIMENSIONS = 3

  type, public :: DistributedMeshForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nDimensions  = 0, &
      nProperCells = 0, &
      nGhostCells  = 0
    integer ( KDI ), private :: &
      iTimer_IO
    integer ( KDI ), dimension ( MAX_N_DIMENSIONS ) :: &
      iaBrick, &
      iaFirst, &
      iaLast, &
      nCells, &
      nGhostLayers, &
      nBricks, &
      nCellsPerBrick
    real ( KDR ) :: &
      CellVolume
    real ( KDR ), dimension ( MAX_N_DIMENSIONS ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      CellWidth, &
      CellArea
    type ( Real_1D_Form ), dimension ( MAX_N_DIMENSIONS ) :: &
      Edge
    type ( MeasuredValueForm ), dimension ( MAX_N_DIMENSIONS ) :: &
      CoordinateUnit
    character ( LDL ) :: &
      BoundaryCondition = 'PERIODIC'
    type ( VariableGroupForm ) :: &
      Position
    type ( VariableGroupForm ), dimension ( : ), allocatable :: &
      ExchangeVariableGroup
    type ( VariableGroup_1D_Form ), allocatable :: &
      VariableGroup
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( PortalHeaderForm ) :: &
      PortalHeaderPrevious, &
      PortalHeaderNext
    type ( MessageIncoming_1D_R_Form ), allocatable :: &
      IncomingPrevious, &
      IncomingNext
    type ( MessageOutgoing_1D_R_Form ), allocatable :: &
      OutgoingPrevious, &
      OutgoingNext
    type ( GridImageStreamForm ) :: &
      GridImageStream
    type ( CurveImageForm ) :: &
      CurveImage
    type ( StructuredGridImageForm ) :: &
      GridImage
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetVariablePointer
    procedure, private, pass :: &
      StartGhostExchangeSingle
    procedure, private, pass :: &
      StartGhostExchangeMultiple
    generic :: &
      StartGhostExchange => StartGhostExchangeSingle, StartGhostExchangeMultiple
    procedure, public, pass :: &
      FinishGhostExchange
    procedure, public, pass :: &
      SetImage
    procedure, public, pass :: &
      Write
    final :: &
      Finalize
  end type DistributedMeshForm

    private :: &
      SetBrickDecomposition, &
      SetGhostExchangePortal, &
      SetCellGeometry, &
      ShowParameters

      private :: &
        BrickIndex

   integer ( KDI ), private, parameter :: &
     TAG_IN_PREV  = 999, &
     TAG_IN_NEXT  = 998, &
     TAG_OUT_PREV = TAG_IN_NEXT, &
     TAG_OUT_NEXT = TAG_IN_PREV

contains


  subroutine Initialize ( DM, C, BoundaryConditionOption )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ), optional :: &
      BoundaryConditionOption

    DM % IGNORABILITY =  CONSOLE % INFO_3
    DM % Communicator => C

    call Show ( 'Initializing a DistributedMesh', DM % IGNORABILITY )
    
    call SetBrickDecomposition ( DM )
    call SetGhostExchangePortal ( DM )
    call SetCellGeometry ( DM )
    
    if ( present ( BoundaryConditionOption ) ) &
      DM % BoundaryCondition = BoundaryConditionOption

    call ShowParameters ( DM )

  end subroutine Initialize


  subroutine SetVariablePointer ( DM, Variable_1D, Variable_3D )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Variable_1D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      Variable_3D

    Variable_3D &
      ( DM % iaFirst ( 1 ) : DM % iaLast ( 1 ), &
        DM % iaFirst ( 2 ) : DM % iaLast ( 2 ), &
        DM % iaFirst ( 3 ) : DM % iaLast ( 3 ) ) &
          => Variable_1D
    
  end subroutine SetVariablePointer


  subroutine StartGhostExchangeSingle ( DM, VG )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM
    class ( VariableGroupForm ), intent ( inout ) :: &
      VG

    allocate ( DM % ExchangeVariableGroup ( 1 ) )
    call DM % ExchangeVariableGroup ( 1 ) % Initialize ( VG )

    call StartGhostExchangeMultiple ( DM, DM % ExchangeVariableGroup )
    
  end subroutine StartGhostExchangeSingle


  subroutine StartGhostExchangeMultiple ( DM, VG )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM
    class ( VariableGroupForm ), dimension ( : ), intent ( inout ) :: &
      VG

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension, etc.
      iG, &  !-- iGroup
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    integer ( KDI ), dimension ( 3 ) :: &
      oSend, &
      nSend
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable

    allocate ( DM % VariableGroup )
    allocate ( DM % IncomingPrevious )
    allocate ( DM % IncomingNext )
    allocate ( DM % OutgoingPrevious )
    allocate ( DM % OutgoingNext )

    associate &
      ( VG_1D  => DM % VariableGroup, &
        PHP => DM % PortalHeaderPrevious, &
        PHN => DM % PortalHeaderNext )

    call VG_1D % Initialize ( VG )

    !-- Post Receives

    call DM % IncomingPrevious % Initialize &
           ( DM % Communicator, spread ( TAG_IN_PREV, 1, PHP % nSources ), &
             PHP % Source, PHP % nChunksFrom * VG_1D % nVariablesTotal )
    call DM % IncomingNext % Initialize &
           ( DM % Communicator, spread ( TAG_IN_NEXT, 1, PHN % nSources ), &
             PHN % Source, PHN % nChunksFrom * VG_1D % nVariablesTotal )

    call DM % IncomingPrevious % Receive ( )
    call DM % IncomingNext % Receive ( )

    !-- Send to Previous

    call DM % OutgoingPrevious % Initialize &
           ( DM % Communicator, spread ( TAG_OUT_PREV, 1, PHP % nTargets ), &
             PHP % Target, PHP % nChunksTo * VG_1D % nVariablesTotal )

    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      oBuffer = 0
      !-- In setting oSend, note Copy command does not inherit lbound
      oSend = DM % nGhostLayers
      nSend ( iD ) = DM % nGhostLayers ( iD )
      nSend ( jD ) = DM % nCellsPerBrick ( jD )
      nSend ( kD ) = DM % nCellsPerBrick ( kD )

      do iG = 1, VG_1D % nGroups
        do iS = 1, VG_1D % nVariables ( iG )          
          iV = VG_1D % VariableGroup ( iG ) % iaSelected ( iS )
          call DM % SetVariablePointer &
                 ( VG_1D % VariableGroup ( iG ) % Value ( :, iV ), V ) 
          call Copy ( V, nSend, oSend, oBuffer, &
                      DM % OutgoingPrevious % Message ( iD ) % Value )
          oBuffer = oBuffer + product ( nSend )
        end do !-- iS
      end do !-- iG

      call DM % OutgoingPrevious % Send ( iD )

    end do !-- iD

    !--- Send to Next

    call DM % OutgoingNext % Initialize &
           ( DM % Communicator, spread ( TAG_OUT_NEXT, 1, PHN % nTargets ), &
             PHN % Target, PHN % nChunksTo * VG_1D % nVariablesTotal )

    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      oBuffer = 0
      !-- In setting oSend, note Copy command does not inherit lbound
      oSend = DM % nGhostLayers
      oSend ( iD ) = oSend ( iD ) + DM % nCellsPerBrick ( iD ) &
                     - DM % nGhostLayers ( iD )
      nSend ( iD ) = DM % nGhostLayers ( iD )
      nSend ( jD ) = DM % nCellsPerBrick ( jD )
      nSend ( kD ) = DM % nCellsPerBrick ( kD )

      do iG = 1, VG_1D % nGroups
        do iS = 1, VG_1D % nVariables ( iG )          
          iV = VG_1D % VariableGroup ( iG ) % iaSelected ( iS )
          call DM % SetVariablePointer &
                 ( VG_1D % VariableGroup ( iG ) % Value ( :, iV ), V ) 
          call Copy ( V, nSend, oSend, oBuffer, &
                      DM % OutgoingNext % Message ( iD ) % Value )
          oBuffer = oBuffer + product ( nSend )
        end do !-- iS
      end do !-- iG

      call DM % OutgoingNext % Send ( iD )

    end do !-- iD

    nullify ( V )
    end associate !-- VG_1D, etc.

  end subroutine StartGhostExchangeMultiple


  subroutine FinishGhostExchange ( DM )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension, etc.
      iG, &  !-- iGroup
      iS, &  !-- iSelected
      iV, &  !-- iVariable
      oBuffer
    integer ( KDI ), dimension ( 3 ) :: &
      oReceive, &
      nReceive
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      V  !-- Variable

    associate ( VG_1D  => DM % VariableGroup )

    !-- Receive from Next

    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      oBuffer = 0
      !-- In setting oReceive, note Copy command does not inherit lbound
      oReceive = DM % nGhostLayers
      oReceive ( iD ) = oReceive ( iD ) + DM % nCellsPerBrick ( iD )
      nReceive ( iD ) = DM % nGhostLayers ( iD )
      nReceive ( jD ) = DM % nCellsPerBrick ( jD )
      nReceive ( kD ) = DM % nCellsPerBrick ( kD )

      call DM % IncomingNext % Wait ( iD )

      do iG = 1, VG_1D % nGroups
        do iS = 1, VG_1D % nVariables ( iG )          
          iV = VG_1D % VariableGroup ( iG ) % iaSelected ( iS )
          call DM % SetVariablePointer &
                 ( VG_1D % VariableGroup ( iG ) % Value ( :, iV ), V ) 
          call Copy ( DM % IncomingNext % Message ( iD ) % Value, &
                      nReceive, oReceive, oBuffer, V )
          oBuffer = oBuffer + product ( nReceive )
        end do !-- iS
      end do !-- iG

    end do !-- iD

    !-- Receive from Previous

    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      oBuffer = 0
      !-- In setting oReceive, note Copy command does not inherit lbound
      oReceive = DM % nGhostLayers
      oReceive ( iD ) = oReceive ( iD ) - DM % nGhostLayers ( iD )
      nReceive ( iD ) = DM % nGhostLayers ( iD )
      nReceive ( jD ) = DM % nCellsPerBrick ( jD )
      nReceive ( kD ) = DM % nCellsPerBrick ( kD )

      call DM % IncomingPrevious % Wait ( iD )

      do iG = 1, VG_1D % nGroups
        do iS = 1, VG_1D % nVariables ( iG )          
          iV = VG_1D % VariableGroup ( iG ) % iaSelected ( iS )
          call DM % SetVariablePointer &
                 ( VG_1D % VariableGroup ( iG ) % Value ( :, iV ), V ) 
          call Copy ( DM % IncomingPrevious % Message ( iD ) % Value, &
                      nReceive, oReceive, oBuffer, V )
          oBuffer = oBuffer + product ( nReceive )
        end do !-- iS
      end do !-- iG

    end do !-- iD

    !-- Cleanup

    call DM % OutgoingPrevious % Wait ( )
    call DM % OutgoingNext % Wait ( )

    nullify ( V )
    end associate !-- VG_1D

    deallocate ( DM % OutgoingNext )
    deallocate ( DM % OutgoingPrevious )
    deallocate ( DM % IncomingNext )
    deallocate ( DM % IncomingPrevious )
    deallocate ( DM % VariableGroup )

    if ( allocated ( DM % ExchangeVariableGroup ) ) &
      deallocate ( DM % ExchangeVariableGroup )

  end subroutine FinishGhostExchange


  subroutine SetImage ( DM, VG, Name )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM
    class ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      VG
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iVG    !-- iVariableGroup
    integer ( KDI ), dimension ( MAX_N_DIMENSIONS ) :: &
      nCells, &
      nGhostInner, &
      nGhostOuter, &
      nExteriorInner, &
      nExteriorOuter
    type ( Real_1D_Form ), dimension ( MAX_N_DIMENSIONS ) :: &
      Edge
    character ( LDF ) :: &
      OutputDirectory

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )
    
    call PROGRAM_HEADER % AddTimer ( 'InputOutput', DM % iTimer_IO )
    
    call PROGRAM_HEADER % Timer ( DM % iTimer_IO ) % Start ( )

    associate ( GIS => DM % GridImageStream )
    call GIS % Initialize &
           ( Name, CommunicatorOption = DM % Communicator, &
             WorkingDirectoryOption = OutputDirectory )
    call GIS % Open ( GIS % ACCESS_SET_GRID )

    nGhostInner = DM % nGhostLayers
    nGhostOuter = DM % nGhostLayers
    nExteriorInner = 0
    nExteriorOuter = 0
    where ( DM % iaBrick == 1 )
      nGhostInner = 0
      nExteriorInner = DM % nGhostLayers
    end where
    where ( DM % iaBrick == DM % nBricks )
      nGhostOuter = 0
      nExteriorOuter = DM % nGhostLayers
    end where

    nCells = DM % iaLast - ( DM % iaFirst - 1 ) &
             - nExteriorInner - nExteriorOuter

    !-- shift lbound of Edge to 1
    do iD = 1, DM % nDimensions
      call Edge ( iD ) % Initialize ( size ( DM % Edge ( iD ) % Value ) )
      Edge ( iD ) % Value = DM % Edge ( iD ) % Value 
    end do !-- iD
    
    select case ( DM % nDimensions )
    case ( 1 ) 

      associate ( CI => DM % CurveImage )
      call CI % Initialize ( GIS )
      call CI % SetGrid &
             ( 'Curves', Edge ( 1 ), DM % nProperCells, &
               oValue = nGhostInner ( 1 ) + nExteriorInner ( 1 ), &
               CoordinateUnitOption = DM % CoordinateUnit ( 1 ) )
      do iVG = 1, size ( VG )
        call CI % AddVariableGroup ( VG ( iVG ) )
      end do
      end associate !-- CI

    case default

      associate ( GI => DM % GridImage )
      call GI % Initialize ( GIS ) 
      call GI % SetGrid &
             ( 'Grid', Edge, nCells, nGhostInner, nGhostOuter, &
               nExteriorInner, nExteriorOuter, DM % nDimensions, &
               DM % nProperCells, DM % nGhostCells, &
               CoordinateUnitOption = DM % CoordinateUnit )
      do iVG = 1, size ( VG )
        call GI % AddVariableGroup ( VG ( iVG ) )
      end do
      end associate !-- GI

    end select

    call GIS % Close ( )

    end associate !-- GIS
    
    call PROGRAM_HEADER % Timer ( DM % iTimer_IO ) % Stop ( )

  end subroutine SetImage


  subroutine Write ( DM, TimeOption, CycleNumberOption )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
          
    call PROGRAM_HEADER % Timer ( DM % iTimer_IO ) % Start ( )
    
    associate ( GIS => DM % GridImageStream )
    call GIS % Open ( GIS % ACCESS_CREATE )

    select case ( DM % nDimensions )
    case ( 1 ) 
      call DM % CurveImage % Write &
             ( TimeOption = TimeOption, CycleNumberOption = CycleNumberOption )
    case default
      call DM % GridImage % Write &
             ( TimeOption = TimeOption, CycleNumberOption = CycleNumberOption )
    end select

    call GIS % Close ( )
    end associate !-- GIS
    
    call PROGRAM_HEADER % Timer ( DM % iTimer_IO ) % Stop ( )
    
  end subroutine Write


  subroutine Finalize ( DM )

    type ( DistributedMeshForm ), intent ( inout ) :: &
      DM

    nullify ( DM % Communicator )

  end subroutine Finalize


  subroutine SetBrickDecomposition ( DM )

    class ( DistributedMeshForm ), intent ( inout ) :: &
      DM

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      nGhostLayers, &
      CommunicatorSizeRoot

    DM % nDimensions = 3
    DM % nCells = 32
    call PROGRAM_HEADER % GetParameter &
           ( DM % nCells, 'nCells', nValuesOption = DM % nDimensions )
    if ( DM % nDimensions < 3 ) DM % nCells ( 3 ) = 1
    if ( DM % nDimensions < 2 ) DM % nCells ( 2 ) = 1

    nGhostLayers = 2
    call PROGRAM_HEADER % GetParameter ( nGhostLayers, 'nGhostLayers' )

    DM % nGhostLayers = 0
    DM % nGhostLayers ( 1 : DM % nDimensions ) = nGhostLayers

    CommunicatorSizeRoot &
      = DM % Communicator % Size ** ( 1.0_KDR / DM % nDimensions ) + 0.5_KDR

    if &
      ( CommunicatorSizeRoot ** DM % nDimensions &
        /= DM % Communicator % Size ) &
    then
      call Show &
             ( 'The square root  ( 2D ) or cube root  ( 3D ) of the number ' &
               // 'of MPI processes must be an integer', &
               CONSOLE % ERROR )
      call Show &
             ( DM % Communicator % Size, 'nProcesses', CONSOLE % ERROR )
      call Show &
             ( DM % nDimensions, 'nDimensions', CONSOLE % ERROR )
      call Show &
             ( DM % Communicator % Size ** ( 1.0_KDR / DM % nDimensions ), &
               'nProcessesRoot', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if
    
    do iD = 1, DM % nDimensions
      if &
        ( DM % nCells ( iD ) > 1 &
         .and. mod ( DM % nCells ( iD ), CommunicatorSizeRoot ) /= 0 ) &
      then
        call Show &
               ( 'The square root  ( 2D ) or cube root  ( 3D ) of the number '&
                 // 'of MPI processes must divide evenly into ' &
                 // 'nCells in each dimension', &
                 CONSOLE % ERROR )
        call Show ( iD, 'iDimension', CONSOLE % ERROR )
        call Show &
               ( DM % nCells ( iD ), 'nTopLevelCells', &
                 CONSOLE % ERROR )
        call Show &
               ( DM % Communicator % Size ** ( 1.0_KDR / DM % nDimensions ), &
                 'nProcessesRoot', CONSOLE % ERROR )
        call Show &
               ( ( 1.0_KDR * DM % nCells ( iD ) ) / CommunicatorSizeRoot, &
                 'nBricks', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
    end do
    
    DM % nCellsPerBrick &
      = max ( DM % nCells / CommunicatorSizeRoot, 1 )
    DM % nBricks = DM % nCells / DM % nCellsPerBrick
    DM % iaBrick = BrickIndex ( DM )  

    DM % iaFirst = 1
    DM % iaLast = 1
    DM % iaFirst ( 1 ) = 1 - DM % nGhostLayers ( 1 )
    DM % iaLast ( 1 ) =  DM % nCellsPerBrick ( 1 ) + DM % nGhostLayers ( 1 )
    if ( DM % nDimensions > 1 ) then
      DM % iaFirst ( 2 ) = 1 - DM % nGhostLayers ( 2 )
      DM % iaLast ( 2 ) =  DM % nCellsPerBrick ( 2 ) + DM % nGhostLayers ( 2 )
    end if
    if ( DM % nDimensions > 2 ) then
      DM % iaFirst ( 3 ) = 1 - DM % nGhostLayers ( 3 )
      DM % iaLast ( 3 ) =  DM % nCellsPerBrick ( 3 ) + DM % nGhostLayers ( 3 )
    end if

    DM % nProperCells = product ( DM % nCellsPerBrick )
    DM % nGhostCells  &
      = product ( DM % iaLast - ( DM % iaFirst - 1 ) ) - DM % nProperCells

  end subroutine SetBrickDecomposition


  subroutine SetGhostExchangePortal ( DM )

    class ( DistributedMeshForm ) , intent ( inout )  :: &
      DM

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iB, jB, kB, &  !-- iBrick, jBrick, kBrick
      iD, jD, kD, &  !-- iDimension, etc.
      nSiblingBricks
    integer ( KDI ), dimension ( MAX_N_DIMENSIONS ) :: &
      iaBrick
    integer ( KDI ), dimension ( : ), allocatable :: &
      SiblingBrickRankPrevious, &
      SiblingBrickRankNext, &
      nCellsSendReceive
    integer ( KDI ), dimension ( :, :, : ), allocatable :: &
      Process

    associate ( nB => DM % nBricks )

    allocate ( Process ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

    iP = 0
    do kB = 1, nB ( 3 )
      do jB = 1, nB ( 2 )
        do iB = 1, nB ( 1 )
          Process ( iB, jB, kB ) = iP
          iP = iP + 1
        end do
      end do
    end do

    nSiblingBricks = DM % nDimensions  !-- faces, either previous or next
    
    allocate ( SiblingBrickRankPrevious ( nSiblingBricks ) )
    allocate ( SiblingBrickRankNext ( nSiblingBricks ) )
    allocate ( nCellsSendReceive ( nSiblingBricks ) )

    !-- Face sibling bricks

    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1

      nCellsSendReceive ( iD ) &
        = DM % nGhostLayers ( iD ) &
          * DM % nCellsPerBrick ( jD ) * DM % nCellsPerBrick ( kD )

      !-- Previous
      iaBrick ( iD ) &
        = mod ( DM % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaBrick ( jD ) = DM % iaBrick ( jD )
      iaBrick ( kD ) = DM % iaBrick ( kD )
      SiblingBrickRankPrevious ( iD ) &
        = Process ( iaBrick ( 1 ), iaBrick ( 2 ), iaBrick ( 3 ) )

      !-- Next
      iaBrick ( iD ) = mod ( DM % iaBrick ( iD ), nB ( iD ) ) + 1
      iaBrick ( jD ) = DM % iaBrick ( jD )
      iaBrick ( kD ) = DM % iaBrick ( kD )
      SiblingBrickRankNext ( iD ) &
        = Process ( iaBrick ( 1 ), iaBrick ( 2 ), iaBrick ( 3 ) )

    end do

    call DM % PortalHeaderPrevious % Initialize &
           ( SiblingBrickRankPrevious, SiblingBrickRankPrevious, &
             nCellsSendReceive, nCellsSendReceive )
    call DM % PortalHeaderNext % Initialize &
           ( SiblingBrickRankNext, SiblingBrickRankNext, &
             nCellsSendReceive, nCellsSendReceive )
    call DM % PortalHeaderPrevious % Show &
           ( 'DistributedMesh PortalHeaderPrevious', DM % IGNORABILITY + 1 )
    call DM % PortalHeaderNext % Show &
           ( 'DistributedMesh PortalHeaderNext', DM % IGNORABILITY + 1 )

    end associate !-- nB

  end subroutine SetGhostExchangePortal


  subroutine SetCellGeometry ( DM )

    class ( DistributedMeshForm ) , intent ( inout )  :: &
      DM

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iD, jD, kD, &  !-- iDimension, etc.
      iE, &  !-- iEdge
      iE_1, iE_2, &  !-- indices of first and last edge values
      oE  !-- oEdge (offset)
    type ( Real_1D_Form ), dimension ( MAX_N_DIMENSIONS ) :: &
      GlobalEdge
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Center

    DM % MinCoordinate  = 0.0_KDR
    DM % MaxCoordinate  = 1.0_KDR
    DM % CoordinateUnit = [ UNIT % IDENTITY, UNIT % IDENTITY, UNIT % IDENTITY ]
    call PROGRAM_HEADER % GetParameter &
           ( DM % MinCoordinate, 'MinCoordinate', &
             InputUnitOption = DM % CoordinateUnit )
    call PROGRAM_HEADER % GetParameter &
           ( DM % MaxCoordinate, 'MaxCoordinate', &
             InputUnitOption = DM % CoordinateUnit )
    if ( DM % nDimensions < 3 ) DM % MaxCoordinate ( 3 ) = 0.0_KDR
    if ( DM % nDimensions < 2 ) DM % MaxCoordinate ( 2 ) = 0.0_KDR

    !-- Cell width and edge values

    DM % CellWidth  = 1.0_KDR

    do iD = 1, DM % nDimensions

      !-- Global edges

      iE_1 = 1 - DM % nGhostLayers ( iD )
      iE_2 = DM % nCells ( iD ) + 1 + DM % nGhostLayers ( iD )
    
      call GlobalEdge ( iD ) % Initialize &
             ( nValues = iE_2 - iE_1 + 1, iLowerBoundOption = iE_1 ) 

      DM % CellWidth ( iD ) &
        = ( DM % MaxCoordinate ( iD ) - DM % MinCoordinate ( iD ) )  &
          / DM % nCells ( iD ) 

      do iE = iE_1, iE_2
        GlobalEdge ( iD ) % Value ( iE )  &
          = DM % MinCoordinate ( iD ) + ( iE - 1 ) * DM % CellWidth ( iD )
      end do 

      !-- Local edges

      iE_1 = 1 - DM % nGhostLayers ( iD )
      iE_2 = DM % nCellsPerBrick ( iD ) + 1 + DM % nGhostLayers ( iD )
    
      call DM % Edge ( iD ) % Initialize &
             ( nValues = iE_2 - iE_1 + 1, iLowerBoundOption = iE_1 ) 

      oE = ( DM % iaBrick ( iD ) - 1 ) * DM % nCellsPerBrick ( iD ) 
      do iE = iE_1, iE_2
        DM % Edge ( iD ) % Value ( iE ) &
          = GlobalEdge ( iD ) % Value ( oE + iE )
      end do 

    end do

    !-- Cell area and volume

    DM % CellArea = 1.0_KDR
    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      DM % CellArea ( iD ) = DM % CellWidth ( jD ) * DM % CellWidth ( kD )
    end do

    DM % CellVolume = product ( DM % CellWidth )

    !-- Cell center position

    call DM % Position % Initialize &
           ( [ DM % nProperCells + DM % nGhostCells, 3 ], &
             NameOption = 'Position', &
             VariableOption = [ 'Center_X                       ', &
                                'Center_Y                       ', &
                                'Center_Z                       ' ], &
             UnitOption = DM % CoordinateUnit, ClearOption = .true. )
    
    do iD = 1, DM % nDimensions
      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1
      call DM % SetVariablePointer ( DM % Position % Value ( :, iD ), Center )
      associate ( Edge => DM % Edge ( iD ) % Value )
      do iC = DM % iaFirst ( iD ), DM % iaLast ( iD )
        select case ( iD )
        case ( 1 )
          Center ( iC, :, : ) = 0.5 * ( Edge ( iC ) + Edge ( iC + 1 ) )
        case ( 2 )
          Center ( :, iC, : ) = 0.5 * ( Edge ( iC ) + Edge ( iC + 1 ) )
        case ( 3 )
          Center ( :, :, iC ) = 0.5 * ( Edge ( iC ) + Edge ( iC + 1 ) )
        end select !-- iD
      end do
      end associate !-- Edge
    end do

    nullify ( Center )

  end subroutine SetCellGeometry


  subroutine ShowParameters ( DM )

    class ( DistributedMeshForm ) , intent ( in )  :: &
      DM

    call Show ( 'DistributedMesh Parameters', DM % IGNORABILITY )
    call Show ( DM % nDimensions, 'nDimensions', DM % IGNORABILITY )
    call Show ( DM % nCells, 'nCells', DM % IGNORABILITY )
    call Show ( DM % nBricks, 'nBricks', DM % IGNORABILITY )
    call Show ( DM % nCellsPerBrick, 'nCellsPerBrick', DM % IGNORABILITY )
    call Show ( DM % nProperCells, 'nProperCells', DM % IGNORABILITY )
    call Show ( DM % nGhostCells, 'nGhostCells', DM % IGNORABILITY )
    call Show ( DM % nGhostLayers, 'nGhostLayers', DM % IGNORABILITY )
    call Show ( DM % MinCoordinate, DM % CoordinateUnit, 'MinCoordinate', &
                DM % IGNORABILITY )
    call Show ( DM % MaxCoordinate, DM % CoordinateUnit, 'MaxCoordinate', &
                DM % IGNORABILITY )
    call Show ( DM % Edge ( 1 ) % Value, 'Edge 1', CONSOLE % INFO_7 )
    if ( DM % nDimensions > 1 ) &
      call Show ( DM % Edge ( 2 ) % Value, 'Edge 2', CONSOLE % INFO_7 )
    if ( DM % nDimensions > 2 ) &
      call Show ( DM % Edge ( 3 ) % Value, 'Edge 3', CONSOLE % INFO_7 )
    call Show ( DM % CellWidth, 'CellWidth', CONSOLE % INFO_7 )
    call Show ( DM % CellArea, 'CellArea', CONSOLE % INFO_7 )
    call Show ( DM % CellVolume, 'CellVolume', CONSOLE % INFO_7 )
    call Show ( DM % BoundaryCondition, 'BoundaryCondition', DM % IGNORABILITY )

  end subroutine ShowParameters


  function BrickIndex ( DM )  result ( BI ) 

    class ( DistributedMeshForm ) , intent ( in )  :: &
      DM
    integer ( KDI ) , dimension ( MAX_N_DIMENSIONS )  :: &
      BI

    if ( DM % nCells ( 3 ) > 1 ) then
      BI ( 3 ) &
        = ( DM % Communicator % Rank &
            / ( DM % nBricks ( 1 ) * DM % nBricks ( 2 ) ) ) &
          + 1
    else
      BI ( 3 ) = 1
    end if

    if ( DM % nCells ( 2 ) > 1 ) then
      BI ( 2 ) &
        = ( mod &
              ( DM % Communicator % Rank, &
                DM % nBricks ( 1 ) * DM % nBricks ( 2 ) )  &
            / DM % nBricks ( 1 ) ) &
          + 1
    else
      BI ( 2 ) = 1
    end if

    BI ( 1 ) &
      = mod &
          ( mod &
              ( DM % Communicator % Rank, &
                DM % nBricks ( 1 ) * DM % nBricks ( 2 ) ), &
              DM % nBricks ( 1 ) ) &
        + 1

  end function BrickIndex


end module DistributedMesh_Form
