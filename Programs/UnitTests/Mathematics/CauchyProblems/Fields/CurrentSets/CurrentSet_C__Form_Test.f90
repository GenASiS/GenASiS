program CurrentSet_C__Form_Test

  !-- CurrentSet_Chart__Form_Test

  use Basics
  use Manifolds
  use Streams
  use CurrentSets

  implicit none

  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    Velocity_U_Unit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( Stream_C_Form ), allocatable :: &
    SC
  type ( CurrentSet_C_Form ), allocatable :: &
    CSC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'CurrentSet_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SC )
  call SC % Initialize ( C, GIS )

  Velocity_U_Unit  =  UNIT % SPEED_MKS

  allocate ( CSC )
  call CSC % Initialize ( C, Velocity_U_Unit )
  call CSC % SetStream ( SC )

  call   C % Show ( )
  call CSC % Show ( )
  call  SC % Show ( )

  deallocate ( CSC )
  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program CurrentSet_C__Form_Test
