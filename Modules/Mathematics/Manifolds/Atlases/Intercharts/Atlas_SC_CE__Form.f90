!-- Atlas_SC_CE represents an Atlas with a single chart with an inner
!   boundary at finite radius.

module Atlas_SC_CE__Form

  !-- Atlas_SingleChart_CentralExcision_Form

  use Basics
  use Charts
  use Atlas_SC__Form

  implicit none
  private

  type, public, extends ( Atlas_SC_Form ) :: Atlas_SC_CE_Form
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, public, pass :: &
      CreateChart
    procedure, public, pass :: &
      CreateChart_CE
  end type Atlas_SC_CE_Form

contains


  subroutine InitializeBasic &
               ( A, Name, CommunicatorOption, IncludeFacesOption, &
                 IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( Atlas_SC_CE_Form ), intent ( inout ) :: &
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
      A % Type = 'an Atlas_SC_CE' 
    end if

    call A % Atlas_SC_Form % Initialize &
           ( Name, CommunicatorOption, IncludeFacesOption, &
             IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
             iDimensionalityOption )

    call A % SetBoundaryConditionsFace &
           ( [ 'OUTFLOW', 'INFLOW ' ], iDimension = 1 )
    if ( A % nDimensions > 1 ) &
      call A % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    if ( A % nDimensions > 2 ) &
      call A % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 3 )

  end subroutine InitializeBasic


  subroutine CreateChart &
               ( A, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, &
                 nDimensionsOption, nEqualOption )

    class ( Atlas_SC_CE_Form ), intent ( inout ) :: &
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

    call Show ( 'For Atlas_SC_CE_Form, please call CreateChart_CE instead.', &
                CONSOLE % ERROR )
    call Show ( 'CreateChart', 'subroutine', CONSOLE % ERROR )
    call Show ( 'Atlas_SC_CE__Form', 'module', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine CreateChart


  subroutine CreateChart_CE &
               ( A, CoordinateUnitOption, RadiusMaxOption, &
                 RadiusExcisionOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption )

    class ( Atlas_SC_CE_Form ), intent ( inout ) :: &
      A
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    allocate ( Chart_SLD_CE_Form :: A % Chart )

    select type ( C => A % Chart )
    class is ( Chart_SLD_CE_Form )
      call C % Initialize &
             ( A, iChart = 1, &
               CoordinateUnitOption = CoordinateUnitOption, &
               RadiusMaxOption = RadiusMaxOption, &
               RadiusExcisionOption = RadiusExcisionOption, &
               RadialRatioOption = RadialRatioOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nCellsPolarOption = nCellsPolarOption )
    end select !-- C

  end subroutine CreateChart_CE


end module Atlas_SC_CE__Form
