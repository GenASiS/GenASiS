program Geometry_F_CH__Form_Test

  !-- Geometry_Flat_ChartHeader__Form_Test

  use Basics
  use Manifolds
  use Streams
  use Geometries

  implicit none

  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( Stream_CH_Form ), allocatable :: &
    SC
  type ( Geometry_F_CH_Form ), allocatable :: &
    GC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F_CH__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  CoordinateUnit  =  UNIT % KILOMETER

  allocate ( C )
  call C % Initialize_H &
         ( PeriodicOption = [ .true., .true., .true. ], &
           CoordinateUnitOption = CoordinateUnit )

  allocate ( SC )
  call SC % Initialize_H ( C, GIS )

  allocate ( GC )
  call GC % Initialize_H ( C )
  call GC % SetStream ( SC )

  call  C % Show ( )
  call GC % Show ( )
  call SC % Show ( )

  deallocate ( GC )
  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_CH__Form_Test
