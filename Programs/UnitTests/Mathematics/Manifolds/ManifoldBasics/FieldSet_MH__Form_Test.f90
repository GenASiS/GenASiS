program FieldSet_MH__Form_Test

  use Basics
  use ManifoldBasics

  implicit none

  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    FSM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_MH__Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( FSM )
  call FSM % Initialize ( M, 'Fields' ) 

  call   M % Show ( )
  call FSM % Show ( )

  deallocate ( FSM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_MH__Form_Test
