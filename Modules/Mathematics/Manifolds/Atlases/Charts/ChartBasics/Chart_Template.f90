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
      Edge, &
      Center, &
      WidthLeft, &
      WidthRight
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
      CoordinateLabel => null ( ), &   
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
      SetCenterRectangular, &
      SetCenterCylindrical_1, &
      SetCenterSpherical_1, &
      SetCenterSpherical_2, &
      SetWidthLeftRight

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

    allocate ( C % Edge       ( MAX_DIMENSIONS ) )
    allocate ( C % Center     ( MAX_DIMENSIONS ) )
    allocate ( C % WidthLeft  ( MAX_DIMENSIONS ) )
    allocate ( C % WidthRight ( MAX_DIMENSIONS ) )

    C % IsDistributed = Atlas % IsDistributed
    if ( present ( IsDistributedOption ) ) &
      C % IsDistributed = IsDistributedOption

    allocate ( C % IsPeriodic ( MAX_DIMENSIONS ) )
    C % IsPeriodic = .false.
    C % IsPeriodic ( : nD ) = IsPeriodic ( : nD )

    allocate ( C % CoordinateSystem )
    C % CoordinateSystem = 'RECTANGULAR'
    if ( present ( CoordinateSystemOption ) ) &
      C % CoordinateSystem = CoordinateSystemOption
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateSystem, 'CoordinateSystem' )

    allocate ( C % CoordinateLabel ( MAX_DIMENSIONS ) )
    C % CoordinateLabel  =  [ 'X', 'Y', 'Z' ]
    select case ( trim ( C % CoordinateSystem ) )
    case ( 'CYLINDRICAL' )
      C % CoordinateLabel  =  [ 'R_Perp', 'Z     ', 'Phi   ' ]
    case ( 'SPHERICAL' )
      C % CoordinateLabel  =  [ 'R    ', 'Theta', 'Phi  ' ]
    end select !-- CoordinateSystem

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
    C % Center           => C_Source % Center
    C % WidthLeft        => C_Source % WidthLeft
    C % WidthRight       => C_Source % WidthRight
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
      deallocate ( C % WidthRight )
      deallocate ( C % WidthLeft )
      deallocate ( C % Center )
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
      nullify ( C % WidthRight )
      nullify ( C % WidthLeft )
      nullify ( C % Center )
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


  subroutine SetGeometryCell ( C, nC, nGL, iD, EdgeValueOption )

    class ( ChartTemplate ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nC, &   !-- nCells
      nGL, &  !-- nGhostLayers
      iD      !-- iDimension
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      EdgeValueOption

    integer ( KDI ) :: &
      iC    !-- iCell
    real ( KDL ) :: &
      Width_IG, &
      Width_OG

    if ( .not. C % AllocatedValues ) &
      return

    if ( .not. allocated ( C % Edge ( iD ) % Value ) ) &
      call C % Edge ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL + 1, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( C % Center ( iD ) % Value ) ) &
      call C % Center ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( C % WidthLeft ( iD ) % Value ) ) &
      call C % WidthLeft ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( C % WidthRight ( iD ) % Value ) ) &
      call C % WidthRight ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )

    !-- Edge, proper cells
    if ( present ( EdgeValueOption ) ) then
      C % Edge ( iD ) % Value ( 1 : nC + 1 )  =  EdgeValueOption
      C % MinCoordinate ( iD )  =  EdgeValueOption ( 1 )
      C % MaxCoordinate ( iD )  =  EdgeValueOption ( nC + 1 )
    else
      select case ( trim ( C % Spacing ( iD ) ) )
      case ( 'EQUAL' )
        call SetEdgeEqual &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), nC )
      case ( 'GEOMETRIC' )
        call SetEdgeGeometric &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), &
                 C % Ratio ( iD ), nC )
      case ( 'COMPACTIFIED' )
        call SetEdgeCompactified &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % Scale ( iD ), nC )
        C % MinCoordinate ( iD )  =  C % Edge ( iD ) % Value ( 1 )
        C % MaxCoordinate ( iD )  =  C % Edge ( iD ) % Value ( nC + 1 )
      case ( 'PROPORTIONAL' )
        call SetEdgeProportional &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % Ratio ( iD ), &
                 C % Scale ( iD ), nC, C % nEqual )
        C % MaxCoordinate ( iD )  =  C % Edge ( iD ) % Value ( nC + 1 )
      case default
        call Show ( 'Spacing not recognized', CONSOLE % ERROR )
        call Show ( 'Chart_Template', 'module', CONSOLE % ERROR )
        call Show ( 'SetGeometryCell', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end if

    !-- Edge, ghost cells
    associate ( Edge => C % Edge ( iD ) % Value )
    do iC = 1, nGL
      Width_IG  =  Edge ( iC + 1 )       -  Edge ( iC )
      Width_OG  =  Edge ( nC - iC + 2 )  -  Edge ( nC - iC + 1 )
      Edge ( 1 - iC )       =  Edge ( 2 - iC )   -  Width_IG
      Edge ( nC + 1 + iC )  =  Edge ( nC + iC )  +  Width_OG
    end do !-- iC
    end associate !-- Edge

    !-- Center
    select case ( trim ( C % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      call SetCenterRectangular &
             ( C % Center ( iD ) % Value, C % Edge ( iD ) % Value )
    case ( 'CYLINDRICAL' )
      select case ( iD )
      case ( 1 )
        call SetCenterCylindrical_1 &
               ( C % Center ( iD ) % Value, C % Edge ( iD ) % Value, &
                 nC, nGL, iaF = 1 - nGL )
      case ( 2, 3 )
        call SetCenterRectangular &
               ( C % Center ( iD ) % Value, &
                 C % Edge ( iD ) % Value )
      end select
    case ( 'SPHERICAL' )
      select case ( iD )
      case ( 1 )
        call SetCenterSpherical_1 &
               ( C % Center ( iD ) % Value, C % Edge ( iD ) % Value, &
                 nC, nGL, iaF = 1 - nGL )
      case ( 2 )
        call SetCenterSpherical_2 &
               ( C % Center ( iD ) % Value, C % Edge ( iD ) % Value, &
                 nC, nGL, iaF = 1 - nGL )
      case ( 3 )
        call SetCenterRectangular &
               ( C % Center ( iD ) % Value, &
                 C % Edge ( iD ) % Value )
      end select
    case default
      call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
      call Show ( 'Chart_Template', 'module', CONSOLE % ERROR )
      call Show ( 'SetGeometryCell', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    !-- Width
    call SetWidthLeftRight &
           ( C % WidthLeft ( iD ) % Value, C % WidthRight ( iD ) % Value, &
             C % Edge ( iD ) % Value, C % Center ( iD ) % Value )

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


  subroutine SetEdgeEqual ( Edge, MinCoordinate, MaxCoordinate, nC )

    !-- Equal cell widths

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate
    integer ( KDI ), intent ( in ) :: &
      nC

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    Edge ( 1 )  =  MinCoordinate
    Width       =  ( MaxCoordinate - MinCoordinate ) / nC

    do iC = 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
    end do

  end subroutine SetEdgeEqual


  subroutine SetEdgeGeometric &
               ( Edge, MinCoordinate, MaxCoordinate, Ratio, nC )

    !-- Each successive cell width is larger by Ratio

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio
    integer ( KDI ), intent ( in ) :: &
      nC

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    Edge ( 1 )  =  MinCoordinate
    Width       =  ( MaxCoordinate - MinCoordinate ) &
                   * ( Ratio - 1.0_KDR ) / ( Ratio ** nC  -  1.0_KDR )

    do iC = 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio * Width
    end do

  end subroutine SetEdgeGeometric


  subroutine SetEdgeCompactified ( Edge, Scale, nC )

    !-- Compactify the domain [ 0, Infinity ] to [ 0, 1 ] via the
    !   transformation Coordinate = Scale * S / ( 1 - S )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nC

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      dS, &
      S, &
      Width

    dS = 1.0_KDR / nC
    S  = 0.5_KDR * dS

    Edge ( 1 )  =  0.0_KDR
    Width       =  Scale  *  dS / ( 1.0_KDR - S ) ** 2 

    do iC = 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      S            =  ( 0.5_KDR + ( iC - 1 ) ) * dS
      Width        =  Scale  *  dS / ( 1.0_KDR - S ) ** 2
    end do

  end subroutine SetEdgeCompactified


  subroutine SetEdgeProportional &
               ( Edge, MinCoordinate, Ratio, Scale, nC, nEqual )

    !-- Width proportional to the inner edge coordinate of the cell

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      Ratio, &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nC, &
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

    do iC = nEqual + 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio  *  Edge ( iC )
    end do

  end subroutine SetEdgeProportional


  subroutine SetCenterRectangular ( Center, Edge )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, size ( Center )
      Center ( iC )  =  0.5_KDR * ( Edge ( iC )  +  Edge ( iC + 1 ) )
    end do

  end subroutine SetCenterRectangular


  subroutine SetCenterCylindrical_1 ( Center, Edge, nC, nGL, iaF )

    real ( KDR ), dimension ( iaF : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( iaF : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nC, &
      nGL, &
      iaF

    integer ( KDI ) :: &
      iC  !-- iCell

    !-- Proper and outer ghost cells
    do iC = 1, nC + nGL
      Center ( iC )  &
        =  2.0_KDR * ( Edge ( iC + 1 ) ** 3  -  Edge ( iC ) ** 3 )  &
           / ( 3.0_KDR * ( Edge ( iC + 1 ) ** 2  -  Edge ( iC ) ** 2 ) )
    end do

    !-- Inner ghost cells
    do iC = 1, nGL
      if ( 0.5_KDR * ( Edge ( 2 - iC )  +  Edge ( 1 - iC ) )  >  0.0_KDR ) &
      then
        Center ( 1 - iC )  &
          =  2.0_KDR * ( Edge ( 2 - iC ) ** 3  -  Edge ( 1 - iC ) ** 3 )  &
             / ( 3.0_KDR * ( Edge ( 2 - iC ) ** 2  -  Edge ( 1 - iC ) ** 2 ) )
      else
        Center ( 1 - iC )  =  - Center ( iC )
      end if
    end do

  end subroutine SetCenterCylindrical_1


  subroutine SetCenterSpherical_1 ( Center, Edge, nC, nGL, iaF )

    real ( KDR ), dimension ( iaF : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( iaF : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nC, &
      nGL, &
      iaF

    integer ( KDI ) :: &
      iC  !-- iCell

    !-- Proper and outer ghost cells
    do iC = 1, nC + nGL
      Center ( iC )  &
        =  3.0_KDR * ( Edge ( iC + 1 ) ** 4  -  Edge ( iC ) ** 4 )  &
           / ( 4.0_KDR * ( Edge ( iC + 1 ) ** 3  -  Edge ( iC ) ** 3 ) )
    end do

    !-- Inner ghost cells
    do iC = 1, nGL
      if ( 0.5_KDR * ( Edge ( 2 - iC )  +  Edge ( 1 - iC ) )  >  0.0_KDR ) &
      then
        Center ( 1 - iC )  &
          =  3.0_KDR * ( Edge ( 2 - iC ) ** 4  -  Edge ( 1 - iC ) ** 4 )  &
             / ( 4.0_KDR * ( Edge ( 2 - iC ) ** 3  -  Edge ( 1 - iC ) ** 3 ) )
      else
        Center ( 1 - iC )  =  - Center ( iC )
      end if
    end do

  end subroutine SetCenterSpherical_1


  subroutine SetCenterSpherical_2 ( Center, Edge, nC, nGL, iaF )

    real ( KDR ), dimension ( iaF : ), intent ( inout ) :: &
      Center
    real ( KDR ), dimension ( iaF : ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nC, &
      nGL, &
      iaF

    integer ( KDI ) :: &
      iC  !-- iCell

    !-- Proper cells
    do iC = 1, nC
      Center ( iC )  &
        =  ( sin ( Edge ( iC + 1 ) )  -  sin ( Edge ( iC ) )  &
             +  Edge ( iC )      *  cos ( Edge ( iC ) )  &
             -  Edge ( iC + 1 )  *  cos ( Edge ( iC + 1 ) ) )  &
           /  ( cos ( Edge ( iC ) )  -  cos ( Edge ( iC + 1 ) ) )
    end do

    !-- Inner ghost cells
    do iC = 1, nGL
      if ( 0.5_KDR * ( Edge ( 2 - iC )  +  Edge ( 1 - iC ) )  >  0.0_KDR ) &
      then
        Center ( 1 - iC )  &
        =  ( sin ( Edge ( 2 - iC ) )  -  sin ( Edge ( 1 - iC ) )  &
             +  Edge ( 1 - iC )  *  cos ( Edge ( 1 - iC ) )  &
             -  Edge ( 2 - iC )  *  cos ( Edge ( 2 - iC ) ) )  &
           /  ( cos ( Edge ( 1 - iC ) )  -  cos ( Edge ( 2 - iC ) ) )
      else
        Center ( 1 - iC )  =  - Center ( iC )
      end if
    end do

    !-- Outer ghost cells
    do iC = 1, nGL
      if ( 0.5_KDR * ( Edge ( nC + iC )  +  Edge ( nC + 1 + iC ) ) &
           <  CONSTANT % PI ) &
      then
        Center ( nC + iC )  &
        =  ( sin ( Edge ( nC + 1 + iC ) )  -  sin ( Edge ( nC + iC ) )  &
             +  Edge ( nC + iC )  *  cos ( Edge ( nC + iC ) )  &
             -  Edge ( nC + 1 + iC )  *  cos ( Edge ( nC + 1 + iC ) ) )  &
           /  (     Edge ( nC + iC )      *  cos ( Edge ( nC + iC ) )  &
                 -  Edge ( nC + 1 + iC )  *  cos ( Edge ( nC + 1 + iC ) ) )
      else
        Center ( nC + iC )  &
          =  2.0_KDR * CONSTANT % PI  -  Center ( nC - ( iC - 1 ) )  
      end if
    end do

  end subroutine SetCenterSpherical_2


  subroutine SetWidthLeftRight ( WidthLeft, WidthRight, Edge, Center )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      WidthLeft, &
      WidthRight
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Edge, &
      Center

    integer ( KDI ) :: &
      iC  !-- iCell

    do iC = 1, size ( Center )
      WidthLeft ( iC )  =  Center ( iC )    -  Edge ( iC )
      WidthRight ( iC )  =  Edge ( iC + 1 )  -  Center ( iC )
    end do

  end subroutine SetWidthLeftRight


end module Chart_Template
