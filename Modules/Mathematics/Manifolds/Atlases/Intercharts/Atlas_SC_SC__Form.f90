!-- Atlas_SC_SC represents an Atlas with a single chart that uses
!   spherical coordinates in 1D, cylindrical coordinates in 2D, and
!   rectangular coordinates in 3D.

module Atlas_SC_SC__Form

  !-- Atlas_SingleChart_SymmetricCurvilinear_Form

  use Basics
  use Charts
  use Atlas_SC__Form

  implicit none
  private

  type, public, extends ( Atlas_SC_Form ) :: Atlas_SC_SC_Form
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, public, pass :: &
      CreateChart
    procedure, public, pass :: &
      CreateChart_SC
  end type Atlas_SC_SC_Form

contains


  subroutine InitializeBasic &
               ( A, Name, CommunicatorOption, IncludeFacesOption, &
                 IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( Atlas_SC_SC_Form ), intent ( inout ) :: &
      A
    character ( * ), intent ( in )  :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    logical ( KDL ), intent ( in ), optional :: &
      IncludeFacesOption, &
      IncludeEdgesOption
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption, &
      nDimensionsOption, &
      iDimensionalityOption

    if ( .not. associated ( A % Type ) ) then
      allocate ( A % Type )
      A % Type = 'an Atlas_SC_SC' 
    end if

    call A % Atlas_SC_Form % Initialize &
           ( Name, CommunicatorOption, IncludeFacesOption, &
             IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
             iDimensionalityOption )

    select case ( A % nDimensions )
    case ( 1 )
      !-- spherical coordinates
      call A % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
    case ( 2 )
      !-- cylindrical coordinates
      call A % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
      call A % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
    case ( 3 )
      !-- rectangular coordinates
      call A % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 1 )
      call A % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
      call A % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 3 )
    end select !-- nDimensions

  end subroutine InitializeBasic


  subroutine CreateChart &
               ( A, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, &
                 nDimensionsOption, nEqualOption )

    class ( Atlas_SC_SC_Form ), intent ( inout ) :: &
      A
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
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    call Show ( 'For Atlas_SC_SC_Form, please call CreateChart_SC instead.', &
                CONSOLE % ERROR )
    call Show ( 'CreateChart', 'subroutine', CONSOLE % ERROR )
    call Show ( 'Atlas_SC_SC__Form', 'module', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine CreateChart


  subroutine CreateChart_SC &
               ( A, CoordinateUnitOption, RadiusMaxOption, &
                 nGhostLayersOption, nCellsRadiusOption )

    class ( Atlas_SC_SC_Form ), intent ( inout ) :: &
      A
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption

    allocate ( Chart_SLD_SC_Form :: A % Chart )

    select type ( C => A % Chart )
    class is ( Chart_SLD_SC_Form )
      call C % Initialize &
             ( A, iChart = 1, &
               CoordinateUnitOption = CoordinateUnitOption, &
               RadiusMaxOption = RadiusMaxOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nCellsRadiusOption = nCellsRadiusOption )
    end select !-- C

  end subroutine CreateChart_SC


end module Atlas_SC_SC__Form
