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
    GeometryName_1 = 'GeometryFlat_Cartesian_1D', &
    GeometryName_2 = 'GeometryFlat_Cartesian_2D', &
    GeometryName_3 = 'GeometryFlat_Cartesian_3D', &
    GeometryName_4 = 'GeometryFlat_Cylindrical_1D', &
    GeometryName_5 = 'GeometryFlat_Cylindrical_2D', &
    GeometryName_6 = 'GeometryFlat_Cylindrical_3D', &
    GeometryName_7 = 'GeometryFlat_Spherical_1D', &
    GeometryName_8 = 'GeometryFlat_Spherical_2D', &
    GeometryName_9 = 'GeometryFlat_Spherical_3D'

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  call TestGeometry ( GeometryName_1, 'CARTESIAN', 1 )
  call TestGeometry ( GeometryName_2, 'CARTESIAN', 2 )
  call TestGeometry ( GeometryName_3, 'CARTESIAN', 3 )

  call TestGeometry ( GeometryName_4, 'CYLINDRICAL', 1 )
  call TestGeometry ( GeometryName_5, 'CYLINDRICAL', 2 )
  call TestGeometry ( GeometryName_6, 'CYLINDRICAL', 3 )

  call TestGeometry ( GeometryName_7, 'SPHERICAL', 1 )
  call TestGeometry ( GeometryName_8, 'SPHERICAL', 2 )
  call TestGeometry ( GeometryName_9, 'SPHERICAL', 3 )

  deallocate ( PROGRAM_HEADER )

end program GeometryFlat_Form_Test


subroutine TestGeometry ( Name, CoordinateSystem, nDimensions )

  use Basics
  use AtlasBasics
  use ProtoChart_Form
  use GeometryFlat_Form

  implicit none

  character ( * ), intent ( in ) :: &
    Name, &
    CoordinateSystem
  integer ( KDI ), intent ( in ) :: &
    nDimensions

  integer ( KDI ) :: &
    i, &
    nCells = 24, &
    nGhostLayers = 2, &
    nEqual = 8
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  type ( AtlasHeaderForm ) :: &
    A
  type ( ProtoChartForm ) :: &
    PC
  type ( GeometryFlatForm ) :: &
    G

  call A % Initialize ( 'Atlas' )
  select case ( trim ( CoordinateSystem ) )
  case ( 'CARTESIAN' )
    call PC % InitializeTemplate &
           ( A, IsPeriodic = [ .false., .false., .false. ], iChart = 1, &
             CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption &
               = [ 0.0_KDR, 5.0_KDR, 10.0_KDR ]  *  UNIT % KILOMETER % Number, &
             MaxCoordinateOption &
               = [ 4.0_KDR, 9.0_KDR, 14.0_KDR ]  *  UNIT % KILOMETER % Number, &
             nDimensionsOption = nDimensions ) 
    CoordinateUnit ( : nDimensions ) = UNIT % KILOMETER
  case ( 'CYLINDRICAL' )
    call PC % InitializeTemplate &
           ( A, IsPeriodic = [ .false., .false., .true. ], iChart = 1, &
             CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption &
               = [ 0.0_KDR, -5.0_KDR, 0.0_KDR ], &
             MaxCoordinateOption &
               = [ 10.0_KDR, 5.0_KDR, 2.0_KDR * CONSTANT % PI ], &
             nDimensionsOption = nDimensions ) 
  case ( 'SPHERICAL' )
    call PC % InitializeTemplate &
           ( A, IsPeriodic = [ .false., .false., .true. ], iChart = 1, &
             SpacingOption = [ 'PROPORTIONAL', 'EQUAL       ', &
                               'EQUAL       ' ], &
             CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption &
               = [ 0.0_KDR, 0.0_KDR, 0.0_KDR ], &
             MaxCoordinateOption &
               = [ 10.0_KDR, CONSTANT % PI, 2.0_KDR * CONSTANT % PI ], &
             RatioOption = [ CONSTANT % PI / 3 * nEqual, 0.0_KDR, 0.0_KDR ], &
             ScaleOption = [ 1.0_KDR, 0.0_KDR, 0.0_KDR ], &
             nDimensionsOption = nDimensions, &
             nEqualOption = nEqual ) 
  end select

  call G % Initialize &
         ( CoordinateSystem, CoordinateUnit, nCells + 2 * nGhostLayers, &
           NameOption = Name )

  call PC % SetGeometryCell ( nCells, nGL = nGhostLayers, iD = 1 )
  G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ) = PC % WidthLeft ( 1 ) % Value
  G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ) = PC % WidthRight ( 1 ) % Value
  G % Value ( :, G % CENTER_U ( 1 ) ) = PC % Center ( 1 ) % Value
  call Show ( PC % Edge ( 1 ) % Value, CoordinateUnit ( 1 ), 'Edge_1' )
  do i = 1 - nGhostLayers, nCells + nGhostLayers
    print *, PC % Edge ( 1 ) % Value ( i ), &
             PC % Center ( 1 ) % Value ( i ), &
             PC % Edge ( 1 ) % Value ( i + 1 )
  end do

  if ( nDimensions > 1 ) then
    call PC % SetGeometryCell ( nCells, nGL = nGhostLayers, iD = 2 )
    G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ) = PC % WidthLeft ( 2 ) % Value
    G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ) = PC % WidthRight ( 2 ) % Value
    G % Value ( :, G % CENTER_U ( 2 ) ) = PC % Center ( 2 ) % Value
    call Show ( PC % Edge ( 2 ) % Value, CoordinateUnit ( 1 ), 'Edge_2' )
    do i = 1 - nGhostLayers, nCells + nGhostLayers
      print *, PC % Edge ( 2 ) % Value ( i ), &
               PC % Center ( 2 ) % Value ( i ), &
               PC % Edge ( 2 ) % Value ( i + 1 )
    end do
  end if

  if ( nDimensions > 2 ) then
    call PC % SetGeometryCell ( nCells, nGL = nGhostLayers, iD = 3 )
    G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ) = PC % WidthLeft ( 3 ) % Value
    G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ) = PC % WidthRight ( 3 ) % Value
    G % Value ( :, G % CENTER_U ( 3 ) ) = PC % Center ( 3 ) % Value
    call Show ( PC % Edge ( 3 ) % Value, CoordinateUnit ( 1 ), 'Edge_3' )
  end if

  call G % SetMetric &
         ( nDimensions, nValues = nCells, oValue = nGhostLayers )

  call Show ( 'Geometry variables' )
  call Show ( G % Name, 'Name' )
  do i = 1, G % nVariables
    call Show ( G % Value ( :, i ), G % Unit ( i ), G % Variable ( i ) )
  end do

end subroutine TestGeometry
