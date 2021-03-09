module Chart_BH__Form

  !-- Chart_BaseHeader_Form

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none
  private

  type, public, extends ( ChartHeaderForm ) :: Chart_BH_Form
  contains
    procedure, private, pass :: &
      InitializeBasic_BH
    generic, public :: &
      Initialize_BH => InitializeBasic_BH
  end type Chart_BH_Form


contains


  subroutine InitializeBasic_BH &
               ( C, M, IsPeriodic, iChart, CommunicatorOption, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, nDimensionsOption, nEqualOption )

    class ( Chart_BH_Form ), intent ( inout ) :: &
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

    ! integer ( KDI ) :: &
    !   iD  !-- iDimension
    ! character ( 2 ) :: &
    !   ChartNumber

    ! C % IGNORABILITY  =   M % IGNORABILITY
    !     C % Manifold  =>  M
    !       C % iChart  =   iChart

    ! C % AllocatedValues = .true.

    ! if ( .not. associated ( C % Type ) ) then
    !   allocate ( C % Type )
    !   C % Type = 'a Chart'
    ! end if

    ! allocate ( C % Name )
    ! write ( ChartNumber, fmt = '(i2.2)' ) iChart
    ! C % Name = 'Chart_' // ChartNumber // '_' // trim ( M % Name ) 

    ! call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    ! call Show ( C % Name, 'Name', C % IGNORABILITY )

    ! if ( present ( nDimensionsOption ) ) then
    !   C % nDimensions  =  nDimensionsOption
    ! else
    !   C % nDimensions  =  M % nDimensions
    ! end if

    ! call SetCoordinateMetadata &
    !        ( C, IsPeriodic, SpacingOption, CoordinateLabelOption, &
    !          CoordinateSystemOption, CoordinateUnitOption, &
    !          MinCoordinateOption, MaxCoordinateOption, RatioOption, &
    !          ScaleOption, nEqualOption )

    ! call SetCells ( C, nCellsOption, nGhostLayersOption )

    ! call SetDecomposition &
    !        ( C, M, CommunicatorOption, nBricksOption, nBricksCompatibleOption )

    ! do iD = 1, C % nDimensions
    !   call SetCoordinateData ( C, iD )
    ! end do !-- iD

  end subroutine InitializeBasic_BH


end module Chart_BH__Form
