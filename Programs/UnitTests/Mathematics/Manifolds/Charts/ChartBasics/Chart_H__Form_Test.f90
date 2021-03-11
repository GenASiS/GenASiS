program Chart_H__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  logical ( KDL ), dimension ( 3 ) :: &
    IsPeriodic
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

  IsPeriodic  =  .true.

  allocate ( Base )
  allocate ( C_Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call C_Base % Initialize_H ( Base, IsPeriodic, iChart = 1 )

  call Base % Show ( )
  call C_Base % Show ( )

  !-- Fiber

  IsPeriodic  =  .false.

  allocate ( Fiber )
  allocate ( C_Fiber )
  call Fiber % Initialize ( 'Fiber', iDimensionalityOption = 2 )
  call C_Fiber % Initialize_H &
         ( Fiber, IsPeriodic, iChart = 1, &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL' )

  call Fiber % Show ( )
  call C_Fiber % Show ( )

  deallocate ( C_Fiber )
  deallocate ( Fiber )
  deallocate ( C_Base )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program Chart_H__Form_Test
