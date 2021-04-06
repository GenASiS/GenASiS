program Chart_H__Form_Test

  !-- Chart_Header_Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( Manifold_H_Form ), allocatable :: &
    Base, &
    Fiber
  type ( Chart_H_Form ), allocatable :: &
    C_Base, &
    C_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Chart_H__Form_Test', DimensionalityOption = '2D_1D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  !-- Base

  Periodic  =  .true.

  allocate ( Base )
  allocate ( C_Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call C_Base % Initialize_H ( Base, 'Global', Periodic )

  !-- Fiber

  Periodic  =  .false.

  allocate ( Fiber )
  allocate ( C_Fiber )
  call Fiber % Initialize ( 'Fiber', iDimensionalityOption = 2 )
  call C_Fiber % Initialize_H &
         ( Fiber, 'Global', Periodic, &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL' )

  !-- Display and cleanup

  call Base % Show ( )
  call Show ( Base % nCharts,    'nCharts',    Base % IGNORABILITY )
  call Show ( Base % nFieldSets, 'nFieldSets', Base % IGNORABILITY )
  call Show ( Base % nStreams,   'nStreams',   Base % IGNORABILITY )

  call C_Base % Show ( )

  call Fiber % Show ( )
  call Show ( Fiber % nCharts,    'nCharts',    Fiber % IGNORABILITY )
  call Show ( Fiber % nFieldSets, 'nFieldSets', Fiber % IGNORABILITY )
  call Show ( Fiber % nStreams,   'nStreams',   Fiber % IGNORABILITY )

  call C_Fiber % Show ( )

  deallocate ( C_Fiber )
  deallocate ( Fiber )
  deallocate ( C_Base )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program Chart_H__Form_Test
