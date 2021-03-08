program Manifold_Form_Test

  use Basics
  use ManifoldBasics

  type ( ManifoldForm ), allocatable :: &
    M
  type ( FieldsHeader_M_Form ), allocatable :: &
    FM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Manifold_Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( FM )
  call FM % Initialize ( M, 'Fields' ) 

  call M % AddFields ( FM )
  call M % Show ( )

  deallocate ( FM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Manifold_Form_Test
