program FieldsHeader_M__Form_Test

  use Basics
  use ManifoldBasics

  type ( ManifoldHeaderForm ), allocatable :: &
    M
  type ( FieldsHeader_M_Form ), allocatable :: &
    FM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldsHeader_M__Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call M % Show ( )

  allocate ( FM )
  call FM % Initialize ( M, 'Fields' ) 

  deallocate ( FM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldsHeader_M__Form_Test
