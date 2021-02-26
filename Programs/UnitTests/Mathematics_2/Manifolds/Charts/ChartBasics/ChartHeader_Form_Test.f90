program ChartHeader_Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  type ( ManifoldHeaderForm ), allocatable :: &
    Base, &
    Fiber
  type ( ChartHeaderForm ), allocatable :: &
    C_Base, &
    C_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'ChartHeader_Form_Test', DimensionalityOption = '2D_1D' )

  allocate ( Base )
  allocate ( C_Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call C_Base % Initialize ( Base, iChart = 1 )
  call Base % Show ( )
  call C_Base % Show ( )

  allocate ( Fiber )
  allocate ( C_Fiber )
  call Fiber % Initialize ( 'Fiber', iDimensionalityOption = 2 )
  call C_Fiber % Initialize ( Fiber, iChart = 1 )
  call Fiber % Show ( )
  call C_Fiber % Show ( )

  deallocate ( C_Fiber )
  deallocate ( Fiber )
  deallocate ( C_Base )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program ChartHeader_Form_Test
