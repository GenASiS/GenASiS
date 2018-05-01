!-- Chart_SLD_SC represents a distributed single-level chart that uses
!   spherical coordinates in 1D, cylindrical coordinates in 2D, and
!   rectangular coordinates in 3D.

module Chart_SLD_SC__Form

  !-- Chart_SingleLevelDistributed_SymmetricCurvilinear_Form

  use Basics
  use AtlasBasics
  use Chart_SLD__Form

  implicit none
  private

  type, public, extends ( Chart_SLD_Form ) :: Chart_SLD_SC_Form
    integer ( KDI ) :: &
      nCellsRadius
    real ( KDR ) :: &
      RadiusMax
  contains
    procedure, private, pass :: &
      Initialize_SC
    generic, public :: &
      Initialize => Initialize_SC
    procedure, private, pass :: &
      ShowHeader
  end type Chart_SLD_SC_Form


contains


  subroutine Initialize_SC &
               ( C, Atlas, iChart, CoordinateUnitOption, RadiusMaxOption, &
                 nGhostLayersOption, nCellsRadiusOption )

    class ( Chart_SLD_SC_Form ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption

    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic
    character ( LDL ) :: &
      CoordinateSystem

    IsPeriodic = .false.

    C % RadiusMax = 1.0_KDR
    if ( present ( RadiusMaxOption ) ) &
      C % RadiusMax = RadiusMaxOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusMax, 'RadiusMax' )

    C % nCellsRadius = 32
    if ( present ( nCellsRadiusOption ) ) &
      C % nCellsRadius = nCellsRadiusOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsRadius, 'nCellsRadius' )

    select case ( Atlas % nDimensions )
    case ( 1 )  !-- Spherical coordinates

      CoordinateSystem = 'SPHERICAL' 

      MinCoordinate = [ 0.0_KDR,   0.0_KDR, 0.0_KDR ]
      MaxCoordinate = [ C % RadiusMax, 0.0_KDR, 0.0_KDR ]

      nCells = [ C % nCellsRadius, 1, 1 ]

    case ( 2 )  !-- Cylindrical coordinates

      CoordinateSystem = 'CYLINDRICAL' 

      MinCoordinate = [ 0.0_KDR,       - C % RadiusMax, 0.0_KDR ]
      MaxCoordinate = [ C % RadiusMax, + C % RadiusMax, 0.0_KDR ]

      nCells = [ C % nCellsRadius, 2 * C % nCellsRadius, 1 ]

    case ( 3 )  !-- Cartesian coordinates
 
      CoordinateSystem = 'CARTESIAN' 

      MinCoordinate  =  - C % RadiusMax
      MaxCoordinate  =  + C % RadiusMax

      nCells  =  2 * C % nCellsRadius
      
    end select !-- nDimensions

    call C % Chart_SLD_Form % Initialize &
           ( Atlas, IsPeriodic, iChart, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayersOption )

  end subroutine Initialize_SC


  subroutine ShowHeader ( C )

    class ( Chart_SLD_SC_Form ), intent ( in ) :: &
      C

    call C % Chart_SLD_Form % Show ( )

    call Show ( 'Chart_SLD_SC parameters' )
    call Show ( C % nCellsRadius, 'nCellsRadius' )
    call Show ( C % RadiusMax, C % CoordinateUnit ( 1 ), 'RadiusMax' )

  end subroutine ShowHeader


end module Chart_SLD_SC__Form

