program Gravitation_G_C__Form_Test

  !-- Gravitation_Galileo_Chart__Form_Test

  use Basics
  use Mathematics
  use Gravitations

  implicit none

  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( Gravitation_G_C_Form ), allocatable :: &
    GC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Gravitation_G_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  CoordinateUnit  =  UNIT % KILOMETER

  allocate ( C )
  call C % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ], &
           CoordinateUnitOption = CoordinateUnit )

  allocate ( SC )
  call SC % Initialize ( C, GIS )

  allocate ( GC )
  call GC % Initialize ( C )
  call GC % SetStream ( SC )

  call  C % Show ( )
  call GC % Show ( )
  call SC % Show ( )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call SC % Write ( )
  call GIS % Close ( )

  deallocate ( GC )
  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Gravitation_G_C__Form_Test
