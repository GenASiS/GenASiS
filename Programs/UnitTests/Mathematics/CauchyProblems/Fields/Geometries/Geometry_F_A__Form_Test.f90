program Geometry_F_A__Form_Test

  !-- Geometry_F_Atlas__Form_Test

  use Basics
  use Manifolds
  use Streams
  use Geometries

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SG_Form ), allocatable :: &
    A
  type ( Stream_ASG_Form ), allocatable :: &
    SA
  type ( Geometry_F_A_Form ), allocatable :: &
    GA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )
  call GA % SetStream ( SA )

  call  A % Show ( )
  call GA % Show ( )
  call SA % Show ( )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call SA % Write ( )
  call GIS % Close ( )

  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_A__Form_Test
