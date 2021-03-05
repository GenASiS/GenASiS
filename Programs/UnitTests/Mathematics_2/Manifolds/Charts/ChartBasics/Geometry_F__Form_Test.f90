program Geometry_F__Form_Test

  use Basics

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F__Form_Test', DimensionalityOption = '1D' )

  associate ( nP => PROGRAM_HEADER % Communicator % Size )
  if ( nP /= 1 ) then
    call Show ( 'This test can only be run with 1 MPI process.', &
                CONSOLE % ERROR )
    call Show ( nP, 'nProcesses', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )
  end if
  end associate !-- nP

  call TestGeometry ( 'RECTANGULAR' )
  call TestGeometry ( 'CYLINDRICAL' )
  call TestGeometry ( 'SPHERICAL' )

  deallocate ( PROGRAM_HEADER )

end program Geometry_F__Form_Test


subroutine TestGeometry ( CoordinateSystem )

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  character ( * ), intent ( in ) :: &
    CoordinateSystem

  integer ( KDI ) :: &
    i
  real ( KDR ) :: &
    MinCoordinate, &
    MaxCoordinate, &
    MinWidth
  type ( MeasuredValueForm ) :: &
    CoordinateUnit
  type ( ManifoldHeaderForm ) :: &
    M
  type ( FieldHeader_M_Form ) :: &
    GM
  type ( ChartHeaderForm ) :: &
    C
  type ( FieldHeader_C_Form ) :: &
    GC
  type ( Geometry_F_Form ) :: &
    G

  call M % Initialize &
         ( 'Manifold_' // trim ( CoordinateSystem ), &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           nDimensionsOption = 1 )
  call M % Show ( )
  call GM % Initialize ( M, 'Geometry' ) 

  MinCoordinate  =   0.0_KDR  *  UNIT % KILOMETER % Number
  MaxCoordinate  =  10.0_KDR  *  UNIT % KILOMETER % Number
       MinWidth  =   0.1_KDR  *  UNIT % KILOMETER % Number

  CoordinateUnit  =  UNIT % KILOMETER

  select case ( trim ( CoordinateSystem ) )
  case ( 'RECTANGULAR' )
    call C % Initialize &
           ( M, IsPeriodic = [ .false., .false., .false. ], iChart = 1, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = [ CoordinateUnit ], &
             MinCoordinateOption = [ MinCoordinate ], &
             MaxCoordinateOption = [ MaxCoordinate ], &
             nDimensionsOption = 1 ) 
  case ( 'CYLINDRICAL' )
    call C % Initialize &
           ( M, IsPeriodic = [ .false., .false., .true. ], iChart = 1, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = [ CoordinateUnit ], &
             MinCoordinateOption = [ MinCoordinate ], &
             MaxCoordinateOption = [ MaxCoordinate ], &
             nDimensionsOption = 1 ) 
  case ( 'SPHERICAL' )
    call C % Initialize &
           ( M, IsPeriodic = [ .false., .false., .true. ], iChart = 1, &
             SpacingOption = [ 'GEOMETRIC' ], &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = [ CoordinateUnit ], &
             MinCoordinateOption = [ MinCoordinate ], &
             MaxCoordinateOption = [ MaxCoordinate ], &
             ScaleOption = [ MinWidth ], &
             nDimensionsOption = 1 )
  end select

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call C % Show ( )

  call GC % Initialize ( GM, C, 'Geometry' ) 

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  call G % Initialize &
         ( GC, nValues = C % nCells ( 1 ) + 2 * C % nGhostLayers ( 1 ), &
           NameOption = 'Geometry_F' )

  !-- Set coordinate fields for 1D and 1 process
  associate &
    (   Edge => C %   Edge ( 1 ) % Value, &
       Width => C %  Width ( 1 ) % Value, &
      Center => C % Center ( 1 ) % Value )
  G % Value ( :, G % EDGE_I_U_1 )  &
    =  Edge ( C % iaFirst ( 1 ) : C % iaLast ( 1 ) )
  G % Value ( :, G % WIDTH_U_1 )  &
    =  Width ( C % iaFirst ( 1 ) : C % iaLast ( 1 ) )
  G % Value ( :, G % CENTER_U_1 )  &
    =  Center ( C % iaFirst ( 1 ) : C % iaLast ( 1 ) )
  end associate !-- Edge, etc.

  call G % ComputeFromCoordinates ( )

  call Show ( 'Geometry variables' )
  call Show ( G % Name, 'Name' )
  do i = 1, G % nVariables
    call Show ( G % Value ( :, i ), G % Unit ( i ), G % Variable ( i ) )
  end do

  call CONSOLE % SetVerbosity ( 'INFO_1' )

end subroutine TestGeometry
