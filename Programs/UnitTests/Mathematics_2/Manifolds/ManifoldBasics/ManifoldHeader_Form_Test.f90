program ManifoldHeader_Form_Test

  use Basics
  use ManifoldBasics

  type ( ManifoldHeaderForm ), allocatable :: &
    Base, &
    Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'ManifoldHeader_Form_Test', DimensionalityOption = '2D_1D' )

  allocate ( Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call Base % Show ( )

  allocate ( Fiber )
  call Fiber % Initialize ( 'Fiber', iDimensionalityOption = 2 )
  call Fiber % Show ( )

  deallocate ( Fiber )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program ManifoldHeader_Form_Test
