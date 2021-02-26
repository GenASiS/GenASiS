program FieldHeader_M__Form_Test

  use Basics
  use ManifoldBasics

  type ( ManifoldHeaderForm ), allocatable :: &
    M
  type ( FieldHeader_M_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldHeader_M__Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call M % Show ( )

  allocate ( F )
  call F % Initialize ( M, 'Field' ) 

  deallocate ( F )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldHeader_M__Form_Test
