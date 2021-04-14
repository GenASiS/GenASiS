program Stream_GS__Form_Test

  !-- Stream_GridStream__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Grid_S_Form ), allocatable :: &
    G
  type ( FieldSet_GS_Form ), allocatable :: &
    FSG
  type ( Stream_GS_Form ), allocatable :: &
    SG

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_GS__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( G )
  call G % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSG )
  call FSG % Initialize ( G )

  allocate ( SG )
  call SG % Initialize ( G, GIS )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SG % AddFieldSet ( FSG )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   G % Show ( )
  call FSG % Show ( )
  call  SG % Show ( )

  deallocate ( SG )
  deallocate ( G )
  deallocate ( FSG )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Stream_GS__Form_Test
