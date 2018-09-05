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
      RadiusScale, &
      RadialRatio, &  !-- nCellsRadial / nCellsPolar
      MinWidth
    type ( CommunicatorForm ), allocatable :: &
      Communicator_2, &
      Communicator_3
  contains
    procedure, public, pass :: &
      InitializeTemplate_C
    procedure, public, pass :: &
      ShowHeaderTemplate_C
    procedure, public, pass :: &
      FinalizeTemplate_C
    procedure, private, pass :: &
      SetCore
  end type Chart_SLD_C_Template

contains


  subroutine InitializeTemplate_C &
               ( C, Atlas, RadiusMin, iChart, CoordinateUnitOption, &
                 RadiusMaxOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nEqualOption )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    real ( KDR ), intent ( in ) :: &
      RadiusMin
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
 
    C % MinWidth  =  C % RadiusScale  *  CONSTANT % PI / nCellsPolar

    nCells = [ nCellsRadial, 1, 1 ]
    if ( Atlas % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( Atlas % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  C % RadiusScale

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


  subroutine ShowHeaderTemplate_C ( C )

    class ( Chart_SLD_C_Template ), intent ( in ) :: &
      C

    call Show ( 'Chart_SLD_C parameters' )
    call Show ( C % nCellsPolar, 'nCellsPolar', C % IGNORABILITY )
    call Show ( C % RadiusScale, C % CoordinateUnit ( 1 ), 'RadiusScale', &
                C % IGNORABILITY )
    call Show ( C % MinWidth, C % CoordinateUnit ( 1 ), 'MinWidth', &
                C % IGNORABILITY )
    call Show ( C % RadialRatio, 'RadialRatio' )
    call Show ( C % RadiusMax, C % CoordinateUnit ( 1 ), &
                'RadiusMax requested', C % IGNORABILITY )
    call Show ( C % MaxCoordinate ( 1 ), C % CoordinateUnit ( 1 ), &
                'RadiusMax actual', C % IGNORABILITY )

  end subroutine ShowHeaderTemplate_C


  impure elemental subroutine FinalizeTemplate_C ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C

    if ( allocated ( C % Communicator_3 ) ) &
      deallocate ( C % Communicator_3 )
    if ( allocated ( C % Communicator_2 ) ) &
      deallocate ( C % Communicator_2 )

    if ( allocated ( C % nSegmentsTo_3 ) ) &
      deallocate ( C % nSegmentsTo_3 )
    if ( allocated ( C % nSegmentsFrom_3 ) ) &
      deallocate ( C % nSegmentsFrom_3 )
    if ( allocated ( C % nSegmentsTo_2 ) ) &
      deallocate ( C % nSegmentsTo_2 )
    if ( allocated ( C % nSegmentsFrom_2 ) ) &
      deallocate ( C % nSegmentsFrom_2 )

  end subroutine FinalizeTemplate_C


  subroutine SetCore ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C

  end subroutine SetCore


end module Chart_SLD_C__Template
