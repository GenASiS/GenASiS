program Geometry_F_C__Form_Test

  !-- Geometry_Flat_Chart__Form_Test

  use Basics
  use Manifolds
  use Streams
  use Geometries

  implicit none

  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Grid_S_Form ), allocatable :: &
    G
  type ( Stream_GS_Form ), allocatable :: &
    SG
  type ( Geometry_F_C_Form ), allocatable :: &
    GC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F_GS__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  CoordinateUnit  =  UNIT % KILOMETER

  allocate ( G )
  call G % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ], &
           CoordinateUnitOption = CoordinateUnit )

  allocate ( SG )
  call SG % Initialize ( G, GIS )

  allocate ( GC )
  call GC % Initialize ( G )
  call GC % SetStream ( SG )

  call  G % Show ( )
  call GC % Show ( )
  call SG % Show ( )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call SG % Write ( )
  call GIS % Close ( )

  deallocate ( GC )
  deallocate ( SG )
  deallocate ( G )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_C__Form_Test
