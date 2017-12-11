!-- Chart contains the most basic functionality characterizing a
!   coordinate chart.

module Chart_Template

  use Basics
  use AtlasBasics

  implicit none
  private

  type, public, abstract :: ChartTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iChart       = 0, &
      nDimensions  = 0, &
      nFields      = 0, &
      nEqual       = 0
    real ( KDR ), dimension ( : ), pointer :: &
      MinCoordinate => null ( ), &
      MaxCoordinate => null ( ), &
      Ratio => null ( ), &
      Scale => null ( )
    class ( Real_1D_Form ), dimension ( : ), pointer :: &
      Edge
    type ( MeasuredValueForm ), dimension ( : ), pointer :: &
      CoordinateUnit => null ( )
    logical ( KDL ) :: &
      IsDistributed, &
      AllocatedValues = .false.
    logical ( KDL ), dimension ( : ), pointer :: &
      IsPeriodic => null ( )
    character ( LDL ), pointer :: &
      CoordinateSystem => null ( )
    character ( LDL ), dimension ( : ), pointer :: &
      Spacing => null ( )
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( )
    type ( ConnectivityForm ), pointer :: &
      Connectivity => null ( )
    class ( AtlasHeaderForm ), pointer :: &
      Atlas => null ( )
  contains 
    procedure, public, pass :: &
      InitializeTemplateBasic
    procedure, public, pass :: &
      InitializeTemplateClone
    generic, public :: &
      InitializeTemplate => InitializeTemplateBasic, InitializeTemplateClone
    procedure, public, pass :: &
      ShowTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, public, pass :: &
      SetBrick
    procedure, public, pass :: &
      SetGeometryCell
  end type ChartTemplate

  ! type, public :: ChartElementForm
  !   class ( ChartTemplate ), allocatable :: &
  !     Element
  ! contains
  !   final :: &
  !     FinalizeElement
  ! end type ChartElementForm

    private :: &
      BrickIndex, &
      SetEdgeEqual, &
      SetEdgeGeometric, &
      SetEdgeCompactified, &
      SetEdgeProportional, &
      SetCenterCartesian, &
      SetCenterCylindrical_1, &
      SetCenterSpherical_1, &
      SetCenterSpherical_2, &
      SetWidthLeftRight, &
      SetGeometryCellEqual, &
      SetGeometryCellGeometric, &
      SetGeometryCellCompactified, &
      SetGeometryCellProportional

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = ATLAS % MAX_DIMENSIONS

contains


  subroutine InitializeTemplateBasic &
               ( C, Atlas, IsPeriodic, iChart, SpacingOption, &
                 CoordinateSystemOption, IsDistributedOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nDimensionsOption, nEqualOption )

    class ( ChartTemplate ), intent ( inout ) :: &
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
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    character ( 2 ) :: &
      ChartNumber

    C % IGNORABILITY = Atlas % IGNORABILITY
    C % iChart = iChart
    C % Atlas => Atlas

    if ( present ( nDimensionsOption ) ) then
      C % nDimensions = nDimensionsOption
    else
      C % nDimensions = Atlas % nDimensions
    end if

    C % AllocatedValues = .true.

    if ( .not. associated ( C % Type ) ) then
      allocate ( C % Type )
      C % Type = 'a Chart'
    end if

    allocate ( C % Name )
    write ( ChartNumber, fmt = '(i2.2)' ) iChart
    C % Name = 'Chart_' // ChartNumber // '_' // trim ( Atlas % Name ) 

    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    associate ( nD => C % nDimensions )

    allocate ( C % Ratio ( MAX_DIMENSIONS ) )
    C % Ratio = 0.0_KDR
    if ( present ( RatioOption ) ) &
      C % Ratio ( : nD ) = RatioOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Ratio ( : nD ), 'Ratio' )

    allocate ( C % Scale ( MAX_DIMENSIONS ) )
    C % Scale = 0.0_KDR
    if ( present ( ScaleOption ) ) &
      C % Scale ( : nD ) = ScaleOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Scale ( : nD ), 'Scale' )

    C % nEqual = 0
    if ( present ( nEqualOption ) ) &
      C % nEqual = nEqualOption

    allocate ( C % CoordinateUnit ( MAX_DIMENSIONS ) )
    C % CoordinateUnit = [ UNIT % IDENTITY, UNIT % IDENTITY, UNIT % IDENTITY ]
    if ( present ( CoordinateUnitOption ) ) &
      C % CoordinateUnit ( : nD ) = CoordinateUnitOption ( : nD )

    allocate ( C % MinCoordinate ( MAX_DIMENSIONS ) )
    C % MinCoordinate = 0.0_KDR
    if ( present ( MinCoordinateOption ) ) &
      C % MinCoordinate ( : nD ) = MinCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % MinCoordinate ( : nD ), 'MinCoordinate', &
             InputUnitOption = C % CoordinateUnit ( : nD ) )

    allocate ( C % MaxCoordinate ( MAX_DIMENSIONS ) )
    C % MaxCoordinate = 0.0_KDR
    C % MaxCoordinate ( : nD ) = 1.0_KDR
    if ( present ( MaxCoordinateOption ) ) &
      C % MaxCoordinate ( : nD ) = MaxCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % MaxCoordinate ( : nD ), 'MaxCoordinate', &
             InputUnitOption = C % CoordinateUnit ( : nD ) )

    allocate ( C % Edge ( MAX_DIMENSIONS ) )

    C % IsDistributed = Atlas % IsDistributed
    if ( present ( IsDistributedOption ) ) &
      C % IsDistributed = IsDistributedOption

    allocate ( C % IsPeriodic ( MAX_DIMENSIONS ) )
    C % IsPeriodic = .false.
    C % IsPeriodic ( : nD ) = IsPeriodic ( : nD )

    allocate ( C % CoordinateSystem )
    C % CoordinateSystem = 'CARTESIAN'
    if ( present ( CoordinateSystemOption ) ) &
      C % CoordinateSystem = CoordinateSystemOption
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateSystem, 'CoordinateSystem' )

    allocate ( C % Spacing ( MAX_DIMENSIONS ) )
    C % Spacing = ''
    C % Spacing ( : nD ) = 'EQUAL'
    if ( present ( SpacingOption ) ) &
      C % Spacing ( : nD ) = SpacingOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Spacing ( : nD ), 'Spacing' )

    associate ( AC => Atlas % Connectivity )
    allocate ( C % Connectivity )
    call C % Connectivity % Initialize &
           ( C % nDimensions, IncludeFacesOption = ( AC % nFaces > 0 ), &
             IncludeEdgesOption = ( AC % nEdges > 0 ) )
    end associate !-- AC

    end associate !-- nD

  end subroutine InitializeTemplateBasic


  subroutine InitializeTemplateClone ( C, C_Source )

    class ( ChartTemplate ), intent ( inout ) :: &
      C
    class ( ChartTemplate ), intent ( in ), target :: &
      C_Source

    C % IGNORABILITY     =  CONSOLE % INFO_7  !-- NOT COPIED!
    C % iChart           =  C_Source % iChart
    C % nDimensions      =  C_Source % nDimensions
    C % nFields          =  0  !-- NOT COPIED!   
    C % MinCoordinate    => C_Source % MinCoordinate
    C % MaxCoordinate    => C_Source % MaxCoordinate
    C % Edge             => C_Source % Edge
    C % Ratio            => C_Source % Ratio
    C % Scale            => C_Source % Scale
    C % nEqual           =  C_Source % nEqual
    C % CoordinateUnit   => C_Source % CoordinateUnit
    C % IsDistributed    =  C_Source % IsDistributed
    C % AllocatedValues  =  .false.  !-- NOT COPIED!
    C % IsPeriodic       => C_Source % IsPeriodic
    C % CoordinateSystem => C_Source % CoordinateSystem
    C % Spacing          => C_Source % Spacing
    C % Type             => C_Source % Type
    C % Name             => C_Source % Name
    C % Connectivity     => C_Source % Connectivity
    C % Atlas            => C_Source % Atlas

  end subroutine InitializeTemplateClone


  subroutine ShowTemplate ( C )

    class ( ChartTemplate ), intent ( in ) :: &
      C

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( C % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )

    call Show ( C % MinCoordinate, C % CoordinateUnit, 'MinCoordinate', &
                C % IGNORABILITY )
    call Show ( C % MaxCoordinate, C % CoordinateUnit, 'MaxCoordinate', &
                C % IGNORABILITY )
    call Show ( C % CoordinateUnit, 'CoordinateUnit', C % IGNORABILITY )

    call Show ( C % Spacing, 'Spacing', C % IGNORABILITY )

    if ( any ( C % Spacing == 'GEOMETRIC' ) &
         .or. any ( C % Spacing == 'PROPORTIONAL' ) ) &
      call Show ( C % Ratio, 'Ratio', C % IGNORABILITY )

    if ( any ( C % Spacing == 'COMPACTIFIED' ) &
         .or. any ( C % Spacing == 'PROPORTIONAL' ) ) &
      call Show ( C % Scale, C % CoordinateUnit, 'Scale', C % IGNORABILITY )

    if ( any ( C % Spacing == 'PROPORTIONAL' ) ) &
      call Show ( C % nEqual, 'nEqual', C % IGNORABILITY )

    call Show ( C % CoordinateSystem, 'CoordinateSystem', C % IGNORABILITY )

    call Show ( C % IsDistributed, 'IsDistributed', C % IGNORABILITY )
    call Show ( C % IsPeriodic, 'IsPeriodic', C % IGNORABILITY )

  end subroutine ShowTemplate


  impure elemental subroutine FinalizeTemplate ( C )

    class ( ChartTemplate ), intent ( inout ) :: &
      C

    if ( .not. associated ( C % Name ) ) return

    call Show ( 'Finalizing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    if ( C % AllocatedValues ) then
      deallocate ( C % Connectivity )
      deallocate ( C % Spacing )
      deallocate ( C % CoordinateSystem )
      deallocate ( C % IsPeriodic )
      deallocate ( C % Edge )
      deallocate ( C % MaxCoordinate )
      deallocate ( C % MinCoordinate )
      deallocate ( C % CoordinateUnit )
      deallocate ( C % Scale )
      deallocate ( C % Ratio )
      deallocate ( C % Name )
      deallocate ( C % Type )
    else
      nullify ( C % Connectivity )
      nullify ( C % Spacing )
      nullify ( C % CoordinateSystem )
      nullify ( C % IsPeriodic )
      nullify ( C % Edge )
      nullify ( C % MaxCoordinate )
      nullify ( C % MinCoordinate )
      nullify ( C % CoordinateUnit )
      nullify ( C % Scale )
      nullify ( C % Ratio )
      nullify ( C % Name )
      nullify ( C % Type )
    end if !-- AllocatedValues

    nullify ( C % Atlas )

  end subroutine FinalizeTemplate


  ! impure elemental subroutine FinalizeElement ( CE )
    
  !   type ( ChartElementForm ), intent ( inout ) :: &
  !     CE

  !   if ( allocated ( CE % Element ) ) deallocate ( CE % Element )

  ! end subroutine FinalizeElement


  subroutine SetBrick &
               ( C, Communicator, nCells, nCellsBrick, nBricks, iaBrick )

    class ( ChartTemplate ), intent ( in ) :: &
      C
    type ( CommunicatorForm ), intent ( in ) :: &
      Communicator
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells
    integer ( KDI ), dimension ( : ), intent ( out ) :: &
      nCellsBrick, &
      nBricks, &
      iaBrick

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      CommunicatorSizeRoot

    CommunicatorSizeRoot &
      = Communicator % Size ** ( 1.0_KDR / C % nDimensions ) + 0.5_KDR

    nBricks = 1
    nBricks ( : C % nDimensions ) = CommunicatorSizeRoot
    call PROGRAM_HEADER % GetParameter &
           ( nBricks ( : C % nDimensions ), 'nBricks' )
    
    if ( product ( nBricks ) /= Communicator % Size ) then
      call Show ( 'The total number of bricks must equal ' &
                  // 'the number of MPI processes', CONSOLE % ERROR )
      call Show ( Communicator % Size, 'nProcesses', CONSOLE % ERROR )
      call Show ( nBricks ( 1 : C % nDimensions ), 'nBricks', CONSOLE % ERROR )
      call Show ( product ( nBricks ), 'product ( nBricks )', CONSOLE % ERROR )
      call Show ( 'Chart_Template', 'module', CONSOLE % ERROR )
      call Show ( 'SetBrick', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Communicator % Synchronize ( )
      call PROGRAM_HEADER % Abort ( )
    end if

    do iD = 1, C % nDimensions
      if ( mod ( nCells ( iD ), nBricks ( iD ) ) /= 0 ) then
        call Show ( 'nBricks in each dimension must divide evenly into ' &
                    // 'nCells in each dimension', CONSOLE % ERROR )
        call Show ( iD, 'iDimension', CONSOLE % ERROR )
        call Show ( nCells ( iD ), 'nCells', CONSOLE % ERROR )
        call Show ( nBricks ( iD ), 'nBricks', CONSOLE % ERROR )
        call Show ( 'Chart_Template', 'module', CONSOLE % ERROR )
        call Show ( 'SetBrick', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Communicator % Synchronize ( )
        call PROGRAM_HEADER % Abort ( )
      end if
    end do
    
    nCellsBrick = nCells / nBricks
    iaBrick = BrickIndex ( nBricks, nCells, Communicator % Rank )

  end subroutine SetBrick


  subroutine SetGeometryCell ( C, Width, Center, nC, nGL, iD, iaF )

    class ( ChartTemplate ), intent ( inout ) :: &
      C
    real ( KDR ), dimension ( iaF : ), intent ( inout ) :: &
      Width, &
      Center
    integer ( KDI ), intent ( in ) :: &
      nC, &   !-- nCells
      nGL, &  !-- nGhostLayers
      iD, &      !-- iDimension
      iaF    !-- iaFirst

    integer ( KDI ) :: &
      iC    !-- iCell

    select case ( trim ( C % Spacing ( iD ) ) )
    case ( 'EQUAL' )
      call SetGeometryCellEqual &
             ( Width  ( 1 : nC ), Center ( 1 : nC ), &
               C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), nC )
    case ( 'GEOMETRIC' )
      call SetGeometryCellGeometric &
             ( Width  ( 1 : nC ), Center ( 1 : nC ), &
               C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), &
               C % Ratio ( iD ), nC )
    case ( 'COMPACTIFIED' )
      call SetGeometryCellCompactified &
             ( Width  ( 1 : nC ), Center ( 1 : nC ), &
               C % Scale ( iD ), nC )
      C % MinCoordinate ( iD ) = Center ( 1 )  - 0.5_KDR * Width ( 1 )
      C % MaxCoordinate ( iD ) = Center ( nC ) + 0.5_KDR * Width ( nC )
    case ( 'PROPORTIONAL' )
      call SetGeometryCellProportional &
             ( Width  ( 1 : nC ), Center ( 1 : nC ), &
               C % MinCoordinate ( iD ), C % Ratio ( iD ), C % Scale ( iD ), &
               nC, C % nEqual )
      C % MaxCoordinate ( iD ) = Center ( nC ) + 0.5_KDR * Width ( nC )
    case default
      call Show ( 'Spacing type not recognized', CONSOLE % ERROR )
      call Show ( 'Chart_Template', 'module', CONSOLE % ERROR )
      call Show ( 'SetGeometryCell', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    do iC = 1, nGL

      Width ( 1 - iC )  = Width ( iC )
      Width ( nC + iC ) = Width ( nC - ( iC - 1 ) )

      Center ( 1 - iC ) &
        = Center ( 1 - ( iC - 1 ) ) &
          - 0.5_KDR * ( Width ( 1 - ( iC - 1 ) ) + Width ( 1 - iC ) ) 
      Center ( nC + iC ) &
        = Center ( nC + ( iC - 1 ) ) &
          + 0.5_KDR * ( Width ( nC + ( iC - 1 ) ) + Width ( nC + iC ) )
 
    end do !-- iC

  end subroutine SetGeometryCell


  function BrickIndex ( nBricks, nCells, MyRank )  result ( BI ) 

    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nBricks, &
      nCells
    integer ( KDI ), intent ( in ) :: &
      MyRank
    integer ( KDI ) , dimension ( ATLAS % MAX_DIMENSIONS )  :: &
      BI

    associate ( nB => nBricks )

    if ( nCells ( 3 ) > 1 ) then
      BI ( 3 ) = ( MyRank / ( nB ( 1 ) * nB ( 2 ) ) ) + 1
    else
      BI ( 3 ) = 1
    end if

    if ( nCells ( 2 ) > 1 ) then
      BI ( 2 ) = ( mod ( MyRank, nB ( 1 ) * nB ( 2 ) ) / nB ( 1 ) ) + 1
    else
      BI ( 2 ) = 1
    end if

    BI ( 1 ) = mod ( mod ( MyRank, nB ( 1 ) * nB ( 2 ) ), nB ( 1 ) ) + 1

    end associate !-- nB

  end function BrickIndex


  subroutine SetEdgeEqual ( Edge, MinCoordinate, MaxCoordinate, nCells )

    !-- Equal cell widths

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    Edge ( 1 )  =  MinCoordinate
    Width       =  ( MaxCoordinate - MinCoordinate ) / nCells

    do iC = 2, nCells + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
    end do

  end subroutine SetEdgeEqual


  subroutine SetEdgeGeometric &
               ( Edge, MinCoordinate, MaxCoordinate, Ratio, nCells )

    !-- Each successive cell width is larger by Ratio

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    Edge ( 1 )  =  MinCoordinate
    Width       =  ( MaxCoordinate - MinCoordinate ) &
                   * ( Ratio - 1.0_KDR ) / ( Ratio ** nCells  -  1.0_KDR )

    do iC = 2, nCells + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio * Width
    end do

  end subroutine SetEdgeGeometric


  subroutine SetEdgeCompactified ( Edge, Scale, nCells )

    !-- Compactify the domain [ 0, Infinity ] to [ 0, 1 ] via the
    !   transformation Coordinate = Scale * S / ( 1 - S )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      dS, &
      S, &
      Width

    dS = 1.0_KDR / nCells
    S  = 0.5_KDR * dS

    Edge ( 1 )  =  0.0_KDR
    Width       =  Scale  *  dS / ( 1.0_KDR - S ) ** 2 

    do iC = 2, nCells + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      S            =  ( 0.5_KDR + ( iC - 1 ) ) * dS
      Width        =  Scale  *  dS / ( 1.0_KDR - S ) ** 2
    end do

  end subroutine SetEdgeCompactified


  subroutine SetEdgeProportional &
               ( Edge, MinCoordinate, Ratio, Scale, nCells, nEqual )

    !-- Width proportional to the inner edge coordinate of the cell

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      Ratio, &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nCells, &
      nEqual

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    if ( nEqual == 0 ) then
      Edge ( 1 )  =  MinCoordinate
    else
      call SetEdgeEqual ( Edge, MinCoordinate, Scale, nEqual )
    end if
    
    Width  =  Ratio  *  Edge ( nEqual + 1 )

    do iC = nEqual + 2, nCells
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio  *  Edge ( iC )
    end do

  end subroutine SetEdgeProportional


  subroutine SetCenterCartesian ( Center, Edge, nCells )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, nCells
      Center ( iC )  =  0.5_KDR * ( Edge ( iC )  +  Edge ( iC + 1 ) )
    end do

  end subroutine SetCenterCartesian


  subroutine SetCenterCylindrical_1 ( Center, Edge, nCells )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, nCells
      Center ( iC )  &
        =  2.0_KDR * ( Edge ( iC + 1 ) ** 3  -  Edge ( iC ) ** 3 )  &
           / ( 3.0_KDR * ( Edge ( iC + 1 ) ** 2  -  Edge ( iC ) ** 2 ) )
    end do

  end subroutine SetCenterCylindrical_1


  subroutine SetCenterSpherical_1 ( Center, Edge, nCells )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, nCells
      Center ( iC )  &
        =  3.0_KDR * ( Edge ( iC + 1 ) ** 4  -  Edge ( iC ) ** 4 )  &
           / ( 4.0_KDR * ( Edge ( iC + 1 ) ** 3  -  Edge ( iC ) ** 3 ) )
    end do

  end subroutine SetCenterSpherical_1


  subroutine SetCenterSpherical_2 ( Center, Edge, nCells )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, nCells
      Center ( iC )  &
        =  ( sin ( Edge ( iC + 1 ) )  -  sin ( Edge ( iC ) )  &
             +  Edge ( iC )      *  cos ( Edge ( iC ) )  &
             -  Edge ( iC + 1 )  *  cos ( Edge ( iC + 1 ) ) )  &
           /  (     Edge ( iC )      *  cos ( Edge ( iC ) )  &
                 -  Edge ( iC + 1 )  *  cos ( Edge ( iC + 1 ) ) )
    end do

  end subroutine SetCenterSpherical_2


  subroutine SetWidthLeftRight ( Width_L, Width_R, Edge, Center, nCells )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Width_L, &
      Width_R
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge, &
      Center
    integer ( KDI ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, nCells
      Width_L ( iC )  =  Center ( iC )    -  Edge ( iC )
      Width_R ( iC )  =  Edge ( iC + 1 )  -  Center ( iC )
    end do

  end subroutine SetWidthLeftRight


  subroutine SetGeometryCellEqual &
               ( Width, Center, MinCoordinate, MaxCoordinate, nCells )

    !-- Equal cell widths

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Width, &
      Center
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      CellWidth

    CellWidth &
      = ( MaxCoordinate - MinCoordinate ) / nCells

    do iC = 1, nCells
      Width  ( iC ) = CellWidth
      Center ( iC ) = MinCoordinate  +  ( 0.5_KDR + ( iC - 1 ) ) * CellWidth
    end do

  end subroutine SetGeometryCellEqual


  subroutine SetGeometryCellGeometric &
               ( Width, Center, MinCoordinate, MaxCoordinate, Ratio, nCells )

    !-- Each successive cell width is larger by Ratio

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Width, &
      Center
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell

    Width ( 1 ) &
      = ( MaxCoordinate - MinCoordinate ) &
        * ( Ratio - 1.0_KDR ) / ( Ratio ** nCells  -  1.0_KDR )

    Center ( 1 ) = MinCoordinate  +  0.5_KDR * Width ( 1 )

    do iC = 2, nCells
      Width  ( iC ) &
        = Ratio * Width ( iC - 1 )
      Center ( iC ) &
        = Center ( iC - 1 )  +  0.5_KDR * ( Width ( iC - 1 ) + Width ( iC ) )
    end do

  end subroutine SetGeometryCellGeometric


  subroutine SetGeometryCellCompactified ( Width, Center, Scale, nCells )

    !-- Compactify the domain [ 0, Infinity ] to [ 0, 1 ] via the
    !   transformation Coordinate = Scale * S / ( 1 - S )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Width, &
      Center
    real ( KDR ), intent ( in ) :: &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nCells

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      dS, &
      S

    dS = 1.0_KDR / nCells
    S  = 0.5_KDR * dS

    Width ( 1 )  = Scale  *  dS / ( 1.0_KDR - S ) ** 2 
    Center ( 1 ) = 0.5_KDR * Width ( 1 )

    do iC = 2, nCells
      S = ( 0.5_KDR + ( iC - 1 ) ) * dS
      Width  ( iC ) &
        = Scale  *  dS / ( 1.0_KDR - S ) ** 2
      Center ( iC ) &
        = Center ( iC - 1 )  +  0.5_KDR * ( Width ( iC - 1 ) + Width ( iC ) )
    end do

  end subroutine SetGeometryCellCompactified


  subroutine SetGeometryCellProportional &
               ( Width, Center, MinCoordinate, Ratio, Scale, nCells, nEqual )

    !-- Width proportional to the inner edge coordinate of the cell

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Width, &
      Center
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      Ratio, &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nCells, &
      nEqual

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      InnerEdge

    if ( nEqual == 0 ) then
      InnerEdge  =  MinCoordinate
    else
      call SetGeometryCellEqual &
             ( Width, Center, MinCoordinate, Scale, nEqual )
      InnerEdge  =  Scale
    end if
    
    Width  ( nEqual + 1 )  =  Ratio * InnerEdge
    Center ( nEqual + 1 )  =  InnerEdge  +  0.5_KDR * Width ( nEqual + 1 ) 

    do iC = nEqual + 2, nCells
      InnerEdge      =  Center ( iC - 1 )  +  0.5_KDR * Width ( iC - 1 )
      Width  ( iC )  =  Ratio * InnerEdge
      Center ( iC )  =  InnerEdge  +  0.5_KDR * Width ( iC )
    end do

  end subroutine SetGeometryCellProportional


end module Chart_Template
