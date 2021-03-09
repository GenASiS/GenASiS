program FieldSet_MH__Form_Test

  use Basics
  use ManifoldBasics

  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    FM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_MH__Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call M % Show ( )

  allocate ( FM )
  call FM % Initialize ( M, 'Fields' ) 

  deallocate ( FM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_MH__Form_Test
