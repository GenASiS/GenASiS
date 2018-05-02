!-- Chart_SLD_CE represents a distributed single-level chart with an inner
!   boundary at finite radius, using spherical coordinates with proportional
!   spacing.

module Chart_SLD_CE__Form

  !-- Chart_SingleLevelDistributed_CentralExcision_Form

  use Basics
  use AtlasBasics
  use Chart_SLD__Form

  implicit none
  private

  type, public, extends ( Chart_SLD_Form ) :: Chart_SLD_CE_Form
    integer ( KDI ) :: &
      nCellsExcision
    real ( KDR ) :: &
      RadiusMin, &
      RadiusMax, &
      RadialRatio  !-- nCellsRadial / nCellsCore
  contains
    procedure, private, pass :: &
      Initialize_CE
    generic, public :: &
      Initialize => Initialize_CE
    procedure, private, pass :: &
      ShowHeader
  end type Chart_SLD_CE_Form

contains


  subroutine Initialize_CE &
               ( C, Atlas, iChart, CoordinateUnitOption, RadiusMaxOption, &
                 RadiusMinOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsExcisionOption )

    class ( Chart_SLD_CE_Form ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusMinOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsExcisionOption

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

    C % RadiusMin = 1.0_KDR
    if ( present ( RadiusMinOption ) ) &
      C % RadiusMin = RadiusMinOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusMin, 'RadiusMin' )

    C % RadiusMax = 10.0_KDR
    if ( present ( RadiusMaxOption ) ) &
      C % RadiusMax = RadiusMaxOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusMax, 'RadiusMax' )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [ C % RadiusMin, 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ C % RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    C % nCellsExcision = 32  !-- Number of "excised" cells with equal spacing
    if ( present ( nCellsExcisionOption ) ) &
      C % nCellsExcision = nCellsExcisionOption
    call PROGRAM_HEADER % GetParameter &
           ( C % nCellsExcision, 'nCellsExcision' )

    C % RadialRatio = 3
    if ( present ( RadialRatioOption ) ) &
      C % RadialRatio = RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    nCellsRadial    = C % RadialRatio  *  C % nCellsExcision 
                      !-- Aim for RadiusMax
    nCellsPolar     = 3  *  C % nCellsExcision  
                      !-- Close to unit aspect ratio
    nCellsAzimuthal = 2  *  nCellsPolar
 
    nCells = [ nCellsRadial, 1, 1 ]
    if ( Atlas % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( Atlas % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  C % RadiusMin

    call C % Chart_SLD_Form % Initialize &
           ( Atlas, IsPeriodic, iChart, SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayersOption )

  end subroutine Initialize_CE


  subroutine ShowHeader ( C )

    class ( Chart_SLD_CE_Form ), intent ( in ) :: &
      C

    call C % Chart_SLD_Form % Show ( )

    call Show ( 'Chart_SLD_CE parameters' )
    call Show ( C % RadiusMin, C % CoordinateUnit ( 1 ), 'RadiusCore' )
    call Show ( C % nCellsExcision, 'nCellsExcision' )
    call Show ( C % RadiusMin / C % nCellsExcision, &
                C % CoordinateUnit ( 1 ), 'CellWidthMin' )
    call Show ( C % RadialRatio, 'RadialRatio' )
    call Show ( C % RadiusMax, C % CoordinateUnit ( 1 ), &
                'RadiusMax requested' )
    call Show ( C % MaxCoordinate ( 1 ), C % CoordinateUnit ( 1 ), &
                'RadiusMax actual' )

  end subroutine ShowHeader


end module Chart_SLD_CE__Form
