program ChartHeader_Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  logical ( KDL ), dimension ( 3 ) :: &
    IsPeriodic
  type ( ManifoldHeaderForm ), allocatable :: &
    Base, &
    Fiber
  type ( ChartHeaderForm ), allocatable :: &
    C_Base, &
    C_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'ChartHeader_Form_Test', DimensionalityOption = '2D_1D' )

  !-- Base

  IsPeriodic  =  .true.

  allocate ( Base )
  allocate ( C_Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call C_Base % Initialize ( Base, IsPeriodic, iChart = 1 )
  call Base % Show ( )
  call C_Base % Show ( )

  !-- Fiber

  IsPeriodic  =  .false.

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( Fiber )
  allocate ( C_Fiber )
  call Fiber % Initialize ( 'Fiber', iDimensionalityOption = 2 )
  call C_Fiber % Initialize &
         ( Fiber, IsPeriodic, iChart = 1, &
           SpacingOption = [ 'GEOMETRIC' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           MinCoordinateOption = [ MinEnergy ], &
           MaxCoordinateOption = [ MaxEnergy ], &
           ScaleOption = [ MinWidthEnergy ] )
  call Fiber % Show ( )
  call C_Fiber % Show ( )

  deallocate ( C_Fiber )
  deallocate ( Fiber )
  deallocate ( C_Base )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program ChartHeader_Form_Test
