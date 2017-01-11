!-- ChartStream_SL handles I/O for single-level charts.

module ChartStream_SL__Form

  !-- ChartStream_SingleLevel_Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use GeometryFlat_CSL__Form

  implicit none
  private

  type, public :: ChartStream_SL_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = ''
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), allocatable :: &
      GridImage
    class ( ChartHeader_SL_Form ), pointer :: &
      Chart => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddField
    procedure, public, pass :: &
      ClearFields
    procedure, public, pass :: &
      Write
    final :: &
      Finalize
  end type ChartStream_SL_Form

  type, public :: ChartStream_SL_ElementForm
    class ( ChartStream_SL_Form ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type ChartStream_SL_ElementForm

    private :: &
      SetEdgeValues

contains


  subroutine Initialize ( CS, Chart, GridImageStream, Name )

    class ( ChartStream_SL_Form ), intent ( inout ) :: &
      CS
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      Chart
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GridImageStream
    character ( * ), intent ( in ) :: &
      Name

    CS % IGNORABILITY = CONSOLE % INFO_3
    CS % Name = Name

    call Show ( 'Initializing a ChartStream_SL', CS % IGNORABILITY )
    call Show ( CS % Name, 'Name', CS % IGNORABILITY )

    CS % GridImageStream => GridImageStream

    select case ( Chart % nDimensions )
    case ( 1 ) 
      allocate ( CS % CurveImage )
      associate ( CI => CS % CurveImage )
      call CI % Initialize ( GridImageStream )
      end associate !-- CI
    case default
      allocate ( CS % GridImage )
      associate ( GI => CS % GridImage )
      call GI % Initialize ( GridImageStream ) 
      end associate !-- GI
    end select !-- nDimensions

    CS % Chart => Chart 

  end subroutine Initialize


  subroutine AddField ( CS, VG )

    class ( ChartStream_SL_Form ), intent ( inout ) :: &
      CS
    class ( VariableGroupForm ), intent ( in ) :: &
      VG

    if ( allocated ( CS % CurveImage ) ) then
      call CS % CurveImage % AddVariableGroup ( VG )
    else if ( allocated ( CS % GridImage ) ) then
      call CS % GridImage % AddVariableGroup ( VG )
    end if

  end subroutine AddField


  subroutine ClearFields ( CS )

    class ( ChartStream_SL_Form ), intent ( inout ) :: &
      CS

    if ( allocated ( CS % CurveImage ) ) then
      call CS % CurveImage % ClearVariableGroups ( )
    else if ( allocated ( CS % GridImage ) ) then
      call CS % GridImage % ClearVariableGroups ( )
    end if

  end subroutine ClearFields


  subroutine Write ( CS, DirectoryOption, TimeOption, CycleNumberOption )

    class ( ChartStream_SL_Form ), intent ( inout ) :: &
      CS
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption

    integer ( KDI ) :: &
      nProperCells, &
      nGhostCells
    integer ( KDI ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      nCells, &
      nGhostInner, &
      nGhostOuter, &
      nExteriorInner, &
      nExteriorOuter
    type ( Real_1D_Form ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
      Edge
    character ( 2 ) :: &
      ChartNumber
    character ( LDF ) :: &
      Directory

    call Show ( 'Writing a ChartStream_SL', CS % IGNORABILITY )
    call Show ( CS % Name, 'Name', CS % IGNORABILITY )

    associate ( C => CS % Chart )

    if ( C % IsDistributed ) then
      nProperCells = C % nProperCells
      nGhostCells  = C % nGhostCells
      nGhostInner = C % nGhostLayers
      nGhostOuter = C % nGhostLayers
      nExteriorInner = 0
      nExteriorOuter = 0
      where ( C % iaBrick == 1 )
        nGhostInner = 0
        nExteriorInner = C % nGhostLayers
      end where
      where ( C % iaBrick == C % nBricks )
        nGhostOuter = 0
        nExteriorOuter = C % nGhostLayers
      end where
    else  ! .not. IsDistributed
      nProperCells = C % nProperCells
      nGhostCells  = 0
      nGhostInner = 0
      nGhostOuter = 0
      nExteriorInner = C % nGhostLayers
      nExteriorOuter = C % nGhostLayers
    end if !-- IsDistributed

    nCells = C % iaLast - ( C % iaFirst - 1 ) - nExteriorInner - nExteriorOuter

    call SetEdgeValues ( CS, Edge )

    write ( ChartNumber, fmt = '(i2.2)' ) C % iChart
    Directory = 'Chart_' // ChartNumber // '/'
    if ( present ( DirectoryOption ) ) &
      Directory = DirectoryOption

    select case ( C % nDimensions )
    case ( 1 ) 
      associate ( CI => CS % CurveImage )
      call CI % SetGrid &
             ( Directory, Edge ( 1 ), nProperCells, &
               oValue = nGhostInner ( 1 ) + nExteriorInner ( 1 ), &
               CoordinateUnitOption = C % CoordinateUnit ( 1 ) )
      call CI % Write &
             ( TimeOption = TimeOption, &
               CycleNumberOption = CycleNumberOption )
      call CI % ClearGrid ( )
      end associate !-- CI
    case default
      associate ( GI => CS % GridImage )
      call GI % SetGrid &
             ( Directory, Edge, nCells, nGhostInner, nGhostOuter, &
               nExteriorInner, nExteriorOuter, C % nDimensions, nProperCells, &
               nGhostCells, CoordinateUnitOption = C % CoordinateUnit )
      call GI % Write &
             ( TimeOption = TimeOption, &
               CycleNumberOption = CycleNumberOption )
      call GI % ClearGrid ( )
      end associate !-- GI
    end select !-- nDimensions

    end associate !-- C

  end subroutine Write


  impure elemental subroutine Finalize ( CS )

    type ( ChartStream_SL_Form ), intent ( inout ) :: &
      CS

    if ( allocated ( CS % GridImage ) ) &
      deallocate ( CS % GridImage )
    if ( allocated ( CS % CurveImage ) ) &
      deallocate ( CS % CurveImage )

    nullify ( CS % GridImageStream )
    nullify ( CS % Chart )

    if ( CS % Name == '' ) return

    call Show ( 'Finalizing a ChartStream_SL', CS % IGNORABILITY )
    call Show ( CS % Name, 'Name', CS % IGNORABILITY )

  end subroutine Finalize


  impure elemental subroutine FinalizeElement ( CSE )
    
    type ( ChartStream_SL_ElementForm ), intent ( inout ) :: &
      CSE

    if ( allocated ( CSE % Element ) ) deallocate ( CSE % Element )

  end subroutine FinalizeElement


  subroutine SetEdgeValues ( CS, Edge )

    class ( ChartStream_SL_Form ), intent ( inout ) :: &
      CS
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Edge

    integer ( KDI ) :: &
      iD,  &  !-- iDimension
      iLB, &  !-- iLowerBound
      iUB     !-- iUpperBound
    real ( KDR ), dimension ( : ), pointer :: &
      Center_1D, &
      Width_1D
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Center_3D, &
      Width_3D

    associate ( C => CS % Chart )

    if ( C % iFieldGeometry == 0 ) then
      call Show ( 'iFieldGeometry not set', CONSOLE % ERROR )
      call Show ( 'ChartStream_SL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetEdgeValues', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Communicator % Synchronize ( )
      call PROGRAM_HEADER % Abort ( )
    end if

    select type ( G_CSL => C % Field ( C % iFieldGeometry ) % Pointer )
    class is ( GeometryFlat_CSL_Form )

    select type ( G => G_CSL % Field )
    class is ( GeometryFlatForm )

    do iD = 1, C % nDimensions

      call C % SetVariablePointer &
             ( G % Value ( :, G % CENTER ( iD ) ), Center_3D )
      call C % SetVariablePointer &
             ( G % Value ( :, G % WIDTH ( iD ) ), Width_3D )

      iLB = lbound ( Center_3D, dim = iD )
      iUB = ubound ( Center_3D, dim = iD )

      select case ( iD )
      case ( 1 )
        Center_1D => Center_3D ( iLB : iUB, 1, 1 )
        Width_1D  => Width_3D  ( iLB : iUB, 1, 1 )
      case ( 2 )
        Center_1D => Center_3D ( 1, iLB : iUB, 1 )
        Width_1D  => Width_3D  ( 1, iLB : iUB, 1 )
      case ( 3 ) 
        Center_1D => Center_3D ( 1, 1, iLB : iUB )
        Width_1D  => Width_3D  ( 1, 1, iLB : iUB )
      end select !-- iD

      call Edge ( iD ) % Initialize ( size ( Center_1D ) + 1 )
      Edge ( iD ) % Value ( 1 )   = Center_1D ( 1 ) - 0.5_KDR * Width_1D ( 1 )
      Edge ( iD ) % Value ( 2 : ) = Center_1D       + 0.5_KDR * Width_1D

    end do !-- iD

    call Show ( Edge ( 1 ) % Value, 'Edge 1', CONSOLE % INFO_7 )
    if ( C % nDimensions > 1 ) &
      call Show ( Edge ( 2 ) % Value, 'Edge 2', CONSOLE % INFO_7 )
    if ( C % nDimensions > 2 ) &
      call Show ( Edge ( 3 ) % Value, 'Edge 3', CONSOLE % INFO_7 )

    end select !-- G
    end select !-- G_CSL
    end associate !-- C

    nullify ( Center_3D )
    nullify ( Width_3D )
    nullify ( Center_1D )
    nullify ( Width_1D )

  end subroutine SetEdgeValues


end module ChartStream_SL__Form
