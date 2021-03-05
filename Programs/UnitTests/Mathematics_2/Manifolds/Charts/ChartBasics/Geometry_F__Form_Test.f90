program Geometry_F__Form_Test

  use Basics

  character ( LDF ) :: &
    GeometryName_1 = 'Geometry_F_Rectangular_1D', &
    GeometryName_2 = 'Geometry_F_Cylindrical_1D', &
    GeometryName_3 = 'Geometry_F_Spherical_1D'

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'Geometry_F__Form_Test' )

  associate ( nP => PROGRAM_HEADER % Communicator % Size )
  if ( nP /= 1 ) then
    call Show ( 'This test can only be run with 1 MPI process.', &
                CONSOLE % ERROR )
    call Show ( nP, 'nProcesses', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )
  end if
  end associate !-- nP

  call TestGeometry ( GeometryName_1, 'RECTANGULAR' )
  call TestGeometry ( GeometryName_2, 'CYLINDRICAL' )

  ! call TestGeometry ( GeometryName_7, 'SPHERICAL', 1 )
  ! call TestGeometry ( GeometryName_8, 'SPHERICAL', 2 )
  ! call TestGeometry ( GeometryName_9, 'SPHERICAL', 3 )

  deallocate ( PROGRAM_HEADER )

end program Geometry_F__Form_Test


subroutine TestGeometry ( Name, CoordinateSystem )

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  character ( * ), intent ( in ) :: &
    Name, &
    CoordinateSystem

  integer ( KDI ) :: &
    i, &
    nCells = 32, &
    nGhostLayers = 2, &
    nEqual = 8
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
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
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call GM % Initialize ( M, 'Geometry' ) 

  CoordinateUnit ( 1 )  =  UNIT % KILOMETER

  select case ( trim ( CoordinateSystem ) )
  case ( 'RECTANGULAR' )
    call C % Initialize &
           ( M, IsPeriodic = [ .false., .false., .false. ], iChart = 1, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnit, &
             MinCoordinateOption &
               = [ 0.0_KDR, 5.0_KDR, 10.0_KDR ] * UNIT % KILOMETER % Number, &
             MaxCoordinateOption &
               = [ 4.0_KDR, 9.0_KDR, 14.0_KDR ] * UNIT % KILOMETER % Number, &
!             nCellsOption = nCells * [ 1, 1, 1 ], &
!             nGhostLayersOption = nGhostLayers * [ 1, 1, 1 ], &
             nDimensionsOption = 1 ) 
  case ( 'CYLINDRICAL' )
    call C % Initialize &
           ( M, IsPeriodic = [ .false., .false., .true. ], iChart = 1, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnit, &
             MinCoordinateOption &
               = [  0.0_KDR  *  UNIT % KILOMETER % Number,  &
                   -5.0_KDR  *  UNIT % KILOMETER % Number, &
                    0.0_KDR  *  CONSTANT % PI ], &
             MaxCoordinateOption &
               = [ 10.0_KDR  *  UNIT % KILOMETER % Number,  &
                    5.0_KDR  *  UNIT % KILOMETER % Number, &
                    2.0_KDR  *  CONSTANT % PI ], &
             nDimensionsOption = 1 ) 
  ! case ( 'SPHERICAL' )
  !   call PC % InitializeTemplate &
  !          ( A, IsPeriodic = [ .false., .false., .true. ], iChart = 1, &
  !            SpacingOption = [ 'PROPORTIONAL', 'EQUAL       ', &
  !                              'EQUAL       ' ], &
  !            CoordinateSystemOption = CoordinateSystem, &
  !            MinCoordinateOption &
  !              = [ 0.0_KDR, 0.0_KDR, 0.0_KDR ], &
  !            MaxCoordinateOption &
  !              = [ 10.0_KDR, CONSTANT % PI, 2.0_KDR * CONSTANT % PI ], &
  !            RatioOption = [ CONSTANT % PI / 3 * nEqual, 0.0_KDR, 0.0_KDR ], &
  !            ScaleOption = [ 1.0_KDR, 0.0_KDR, 0.0_KDR ], &
  !            nDimensionsOption = nDimensions, &
  !            nEqualOption = nEqual ) 
  end select

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call C % Show ( )

  call GC % Initialize ( GM, C, 'Geometry' ) 

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  call G % Initialize &
         ( GC, nCells + 2 * nGhostLayers, NameOption = Name )

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
