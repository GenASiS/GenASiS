program Geometry_N__Form_Test

  use Basics

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Geometry_N__Form_Test', &
    GeometryName_1 = 'Geometry_N_1D', &
    GeometryName_2 = 'Geometry_N_2D', &
    GeometryName_3 = 'Geometry_N_3D'

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  call TestGeometry ( GeometryName_1, 1 )
  call TestGeometry ( GeometryName_2, 2 )
  call TestGeometry ( GeometryName_3, 3 )

  deallocate ( PROGRAM_HEADER )

end program Geometry_N__Form_Test


subroutine TestGeometry ( Name, nDimensions )

  use Basics
  use Mathematics
  use Geometry_N__Form

  implicit none

  character ( LDF ), intent ( in ) :: &
    Name
  integer ( KDI ), intent ( in ) :: &
    nDimensions

  integer ( KDI ) :: &
    i, &
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
  type ( Chart_SLL_Form ) :: &
    C
  type ( Geometry_N_Form ) :: &
    G

  call A % Initialize ( 'Atlas' )
  call C % Initialize &
         ( A, IsPeriodic = [ .false., .false., .false. ], iChart = 1, &
           MinCoordinateOption &
             = [ 0.0_KDR, 5.0_KDR, 10.0_KDR ]  *  UNIT % KILOMETER % Number, &
           MaxCoordinateOption &
             = [ 4.0_KDR, 9.0_KDR, 14.0_KDR ]  *  UNIT % KILOMETER % Number, &
           nDimensionsOption = nDimensions ) 

  CoordinateSystem = 'CARTESIAN'
  CoordinateUnit ( : nDimensions ) = UNIT % KILOMETER

  call G % Initialize &
         ( CoordinateSystem, CoordinateUnit, nValues, NameOption = Name )

  allocate ( Width ( nValues ) )
  allocate ( Center ( nValues ) )

  call C % SetGeometryCell &
         ( Width, Center, nValues, nGL = 0, iD = 1, iaF = 1 )
  G % Value ( :, G % WIDTH ( 1 ) )  = Width
  G % Value ( :, G % CENTER ( 1 ) ) = Center

  if ( nDimensions > 1 ) then
    call C % SetGeometryCell &
           ( Width, Center, nValues, nGL = 0, iD = 2, iaF = 1 )
    G % Value ( :, G % WIDTH ( 2 ) )  = Width
    G % Value ( :, G % CENTER ( 2 ) ) = Center
  end if

  if ( nDimensions > 2 ) then
    call C % SetGeometryCell &
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
