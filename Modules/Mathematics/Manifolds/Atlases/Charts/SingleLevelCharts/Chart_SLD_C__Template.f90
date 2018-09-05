!-- Chart_SLD_C represents a distributed single-level chart using spherical 
!   coordinates and proportional radial spacing.

module Chart_SLD_C__Template

  !-- Chart_SingleLevelDistributed_Central_Template

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLD__Form

  implicit none
  private

  type, public, extends ( Chart_SLD_Form ), abstract :: Chart_SLD_C_Template
    integer ( KDI ) :: &
      nCellsPolar, &
      nPillars_2, &
      nPillars_3
    integer ( KDI ), dimension ( : ), allocatable :: &
      nSegmentsFrom_2, nSegmentsTo_2, &
      nSegmentsFrom_3, nSegmentsTo_3
    real ( KDR ) :: &
      RadiusMax, &
      RadialRatio, &  !-- nCellsRadial / nCellsPolar
      MinWidth
  contains
    procedure, public, pass :: &
      InitializeTemplate_C
    procedure, private, pass :: &
      SetCore
  end type Chart_SLD_C_Template

contains


  subroutine InitializeTemplate_C &
               ( C, Atlas, RadiusMin, RadiusScale, iChart, &
                 CoordinateUnitOption, RadiusMaxOption, RadialRatioOption, &
                 nGhostLayersOption, nCellsPolarOption, nEqualOption )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    real ( KDR ), intent ( in ) :: &
      RadiusMin, &
      RadiusScale
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption

    integer ( KDI ) :: &
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    CoordinateSystem = 'SPHERICAL'

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    IsPeriodic = .false.
    IsPeriodic ( 3 ) = .true.

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [ RadiusMin,     0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ C % RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    C % nCellsPolar = 128
    if ( present ( nCellsPolarOption ) ) &
      C % nCellsPolar = nCellsPolarOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsPolar, 'nCellsPolar' )

    call C % SetCore ( )

    C % RadialRatio = 1
    if ( present ( RadialRatioOption ) ) &
      C % RadialRatio = RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    nCellsRadial     =  C % RadialRatio * C % nCellsPolar !-- Aim for RadiusMax
    nCellsPolar      =  C % nCellsPolar
    nCellsAzimuthal  =  2 * nCellsPolar
 
    nCells = [ nCellsRadial, 1, 1 ]
    if ( Atlas % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( Atlas % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  RadiusScale

    call C % Chart_SLD_Form % Initialize &
           ( Atlas, IsPeriodic, iChart, SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayersOption, &
             nEqualOption = nEqualOption )

  end subroutine InitializeTemplate_C


  subroutine SetCore ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C

  end subroutine SetCore


end module Chart_SLD_C__Template
