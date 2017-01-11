module ProtoChart_Form

  use Chart_Template

  implicit none
  private

  type, public, extends ( ChartTemplate ) :: ProtoChartForm
  end type ProtoChartForm

end module ProtoChart_Form


program GeometryFlat_Form_Test

  use Basics

  implicit none

  character ( LDF ) :: &
    ProgramName = 'GeometryFlat_Form_Test', &
    GeometryName_1 = 'GeometryFlat_1D', &
    GeometryName_2 = 'GeometryFlat_2D', &
    GeometryName_3 = 'GeometryFlat_3D'

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  call TestGeometry ( GeometryName_1, 1 )
  call TestGeometry ( GeometryName_2, 2 )
  call TestGeometry ( GeometryName_3, 3 )

  deallocate ( PROGRAM_HEADER )

end program GeometryFlat_Form_Test


subroutine TestGeometry ( Name, nDimensions )

  use Basics
  use AtlasBasics
  use ProtoChart_Form
  use GeometryFlat_Form

  implicit none

  character ( LDF ), intent ( in ) :: &
    Name
  integer ( KDI ), intent ( in ) :: &
    nDimensions

  integer ( KDI ) :: &
    i, &
    iD, &  !-- iDimension
    nValues = 4
  real ( KDR ), dimension ( : ), allocatable :: &
    Width, &
    Center
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  character ( LDL ) :: &
    CoordinateSystem
  type ( AtlasHeaderForm ) :: &
    A
  type ( ProtoChartForm ) :: &
    PC
  type ( GeometryFlatForm ) :: &
    G

  call A % Initialize ( 'Atlas' )
  call PC % InitializeTemplate &
         ( A, IsPeriodic = [ .false., .false., .false. ], iChart = 1, &
           MinCoordinateOption = [ 0.0_KDR, 5.0_KDR, 10.0_KDR ], &
           MaxCoordinateOption = [ 4.0_KDR, 9.0_KDR, 14.0_KDR ], &
           nDimensionsOption = nDimensions ) 

  CoordinateSystem = 'CARTESIAN'
  CoordinateUnit ( : nDimensions ) = UNIT % KILOMETER

  call G % Initialize &
         ( CoordinateSystem, CoordinateUnit, nValues, NameOption = Name )

  allocate ( Width ( nValues ) )
  allocate ( Center ( nValues ) )

  call PC % SetGeometryCell &
         ( Width, Center, nValues, nGL = 0, iD = 1, iaF = 1 )
  G % Value ( :, G % WIDTH ( 1 ) )  = Width
  G % Value ( :, G % CENTER ( 1 ) ) = Center

  if ( nDimensions > 1 ) then
    call PC % SetGeometryCell &
           ( Width, Center, nValues, nGL = 0, iD = 2, iaF = 1 )
    G % Value ( :, G % WIDTH ( 2 ) )  = Width
    G % Value ( :, G % CENTER ( 2 ) ) = Center
  end if

  if ( nDimensions > 2 ) then
    call PC % SetGeometryCell &
           ( Width, Center, nValues, nGL = 0, iD = 3, iaF = 1 )
    G % Value ( :, G % WIDTH ( 3 ) )  = Width
    G % Value ( :, G % CENTER ( 3 ) ) = Center
  end if

  call G % SetMetric ( nDimensions, nValues, 0 )

  call Show ( 'Geometry variables' )
  call Show ( G % Name, 'Name' )
  do i = 1, G % nVariables
    call Show ( G % Value ( :, i ), G % Unit ( i ), G % Variable ( i ) )
  end do

end subroutine TestGeometry
