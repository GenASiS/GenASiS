!-- Chart_SLL represents a local (not distributed) single-level chart.

module Chart_SLL__Form

  !-- Chart_SingleLevelLocal_Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SL__Template

  implicit none
  private

  type, public, extends ( Chart_SL_Template ) :: Chart_SLL_Form
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, public, pass ( C ) :: &
      CopyBoundary
    procedure, public, pass ( C ) :: &
      ReverseBoundary
    final :: &
      Finalize
    procedure, public, pass :: &
      SetGeometryWidthCenter
  end type Chart_SLL_Form

contains


  subroutine InitializeBasic &
               ( C, Atlas, IsPeriodic, iChart, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 IsDistributedOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, &
                 nDimensionsOption, nEqualOption )

    class ( Chart_SLL_Form ), intent ( inout ) :: &
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
      C % Type = 'a Chart_SLL'
    end if

    call C % ChartHeader_SL_Form % Initialize &
           ( Atlas, IsPeriodic, iChart, SpacingOption, CoordinateLabelOption, &
             CoordinateSystemOption, IsDistributedOption, &
             CoordinateUnitOption, MinCoordinateOption, MaxCoordinateOption, &
             RatioOption, ScaleOption, nCellsOption, nGhostLayersOption, &
             nDimensionsOption, nEqualOption )

  end subroutine InitializeBasic


  subroutine CopyBoundary ( Value, C, iDimension, iConnection )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    class ( Chart_SLL_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    call C % CopyBoundaryTemplate &
           ( Value, C % nCells, iDimension, iConnection )

  end subroutine CopyBoundary


  subroutine ReverseBoundary ( Value, C, iDimension, iConnection )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    class ( Chart_SLL_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    call C % ReverseBoundaryTemplate &
           ( Value, C % nCells, iDimension, iConnection )

  end subroutine ReverseBoundary


  impure elemental subroutine Finalize ( C )

    type ( Chart_SLL_Form ), intent ( inout ) :: &
      C

    call C % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetGeometryWidthCenter &
               ( C, Center, Width_L, Width_R, iD, EdgeValueOption )

    class ( Chart_SLL_Form ), intent ( inout ) :: &
      C
    type ( Real_1D_Form ), intent ( inout ) :: &
      Center, &
      Width_L, &
      Width_R
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      EdgeValueOption

    associate &
      ( nC  => C % nCells ( iD ), &
        nGL => C % nGhostLayers ( iD ), &
        iaF => C % iaFirst ( iD ), &
        iaL => C % iaLast ( iD ) )

    call Center % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    call Width_L % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )
    call Width_R % Initialize &
           ( iaL - ( iaF - 1 ), iLowerBoundOption = iaF )

    call C % SetGeometryCell ( nC, nGL, iD, EdgeValueOption )

    call Copy ( C % Center ( iD ) % Value, Center % Value )
    call Copy ( C % WidthLeft ( iD ) % Value, Width_L % Value )
    call Copy ( C % WidthRight ( iD ) % Value, Width_R % Value )

    end associate !-- nC, etc.

  end subroutine SetGeometryWidthCenter


end module Chart_SLL__Form
