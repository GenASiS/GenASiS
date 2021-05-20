program Geometry_G__Form_Test

  use Basics

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Geometry_G__Form_Test', &
    GeometryName_1 = 'Geometry_G_1D', &
    GeometryName_2 = 'Geometry_G_2D', &
    GeometryName_3 = 'Geometry_G_3D'

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )

  call CONSOLE % SetVerbosity ( 'INFO_4' )

  call TestGeometry ( GeometryName_1, 1 )
  call TestGeometry ( GeometryName_2, 2 )
  call TestGeometry ( GeometryName_3, 3 )

  deallocate ( PROGRAM_HEADER )

end program Geometry_G__Form_Test


subroutine TestGeometry ( Name, nDimensions )

  use Basics
  use Mathematics
  use Geometry_G__Form

  implicit none

  character ( LDF ), intent ( in ) :: &
    Name
  integer ( KDI ), intent ( in ) :: &
    nDimensions

  integer ( KDI ) :: &
    i, &
    nValues = 4
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  character ( LDL ) :: &
    CoordinateSystem
  type ( AtlasHeaderForm ) :: &
    A
  type ( Chart_SLL_Form ) :: &
    C
  type ( Geometry_G_Form ) :: &
    G

  call A % Initialize ( 'Atlas' )
  call C % Initialize &
         ( A, IsPeriodic = [ .false., .false., .false. ], iChart = 1, &
           MinCoordinateOption &
             = [ 0.0_KDR, 5.0_KDR, 10.0_KDR ]  *  UNIT % KILOMETER % Number, &
           MaxCoordinateOption &
             = [ 4.0_KDR, 9.0_KDR, 14.0_KDR ]  *  UNIT % KILOMETER % Number, &
           nDimensionsOption = nDimensions ) 

  CoordinateSystem = 'RECTANGULAR'
  CoordinateUnit ( : nDimensions ) = UNIT % KILOMETER

  call G % Initialize &
         ( CoordinateSystem, CoordinateUnit, nValues, NameOption = Name )

  call C % SetGeometryCell ( nValues, nGL = 0, iD = 1 )
  G % Value ( :, G % WIDTH_LEFT_U ( 1 ) )  =  C % WidthLeft ( 1 ) % Value
  G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) )  =  C % WidthRight ( 1 ) % Value
  G % Value ( :, G % CENTER_U ( 1 ) )  =  C % Center ( 1 ) % Value

  if ( nDimensions > 1 ) then
    call C % SetGeometryCell ( nValues, nGL = 0, iD = 2 )
    G % Value ( :, G % WIDTH_LEFT_U ( 2 ) )  =  C % WidthLeft ( 2 ) % Value
    G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) )  =  C % WidthRight ( 2 ) % Value
    G % Value ( :, G % CENTER_U ( 2 ) )  =  C % Center ( 2 ) % Value
  end if

  if ( nDimensions > 2 ) then
    call C % SetGeometryCell ( nValues, nGL = 0, iD = 3 )
    G % Value ( :, G % WIDTH_LEFT_U ( 3 ) )  =  C % WidthLeft ( 3 ) % Value
    G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) )  =  C % WidthRight ( 3 ) % Value
    G % Value ( :, G % CENTER_U ( 3 ) ) =  C % Center ( 3 ) % Value
  end if

  call G % SetMetricFixed ( nDimensions, nValues, 0 )

  call Show ( 'Geometry variables' )
  call Show ( G % Name, 'Name' )
  do i = 1, G % nVariables
    call Show ( G % Value ( :, i ), G % Unit ( i ), G % Variable ( i ) )
  end do

end subroutine TestGeometry
