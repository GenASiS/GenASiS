!-- Chart_Form is a skeletal representation of a Manifold.

module Chart_Form

  use Basics
  use ManifoldBasics
  use ChartHeader_Form
  use FieldsHeader_C__Form

  implicit none
  private

  type, public, extends ( ChartHeaderForm ) :: ChartForm
    integer ( KDI ) :: &
      nFieldSets = 0
    type ( FieldsHeader_C_Pointer ), dimension ( : ), allocatable :: &
      Fields
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, public, pass :: &
      AddFields
    procedure, private, pass :: &
      Show_C
    final :: &
      Finalize
  end type ChartForm

    integer ( KDI ), private, parameter :: &
      MAX_FIELDS  =  MANIFOLD % MAX_FIELDS

contains


  subroutine InitializeBasic &
               ( C, M, IsPeriodic, iChart, CommunicatorOption, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, nDimensionsOption, nEqualOption )

    class ( ChartForm ), intent ( inout ) :: &
      C
    class ( ManifoldHeaderForm ), intent ( in ), target :: &
      M
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsPeriodic
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption, &
      nBricksOption, &
      nBricksCompatibleOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    call C % ChartHeaderForm % Initialize &
           ( M, IsPeriodic, iChart, CommunicatorOption, SpacingOption, &
             CoordinateLabelOption, CoordinateSystemOption, &
             CoordinateUnitOption, MinCoordinateOption, &
             MaxCoordinateOption, RatioOption, ScaleOption, &
             nCellsOption, nGhostLayersOption, nBricksOption, &
             nBricksCompatibleOption, nDimensionsOption, nEqualOption )

    allocate ( C % Fields ( MAX_FIELDS ) )

  end subroutine InitializeBasic


  subroutine AddFields ( C, Fields )

    class ( ChartForm ), intent ( inout ) :: &
      C
    class ( FieldsHeader_C_Form ), intent ( in ), target :: &
      Fields

    associate ( nF  =>  C % nFieldSets )

    nF  =  nF + 1

    C % Fields ( nF ) % Pointer  =>  Fields

    end associate !-- nF

  end subroutine AddFields


  subroutine Show_C ( C )

    class ( ChartForm ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iF

    call C % ChartHeaderForm % Show ( )

    call Show ( C % nFieldSets, 'nFieldSets', C % IGNORABILITY )
    call Show ( [ ( C % Fields ( iF ) % Pointer % Name, &
                    iF = 1, C % nFieldSets) ], &
                'FieldSets', C % IGNORABILITY )

  end subroutine Show_C


  subroutine Finalize ( C )

    type ( ChartForm ), intent ( inout ) :: &
      C

    if ( allocated ( C % Fields ) ) &
      deallocate ( C % Fields )

  end subroutine Finalize


end module Chart_Form
