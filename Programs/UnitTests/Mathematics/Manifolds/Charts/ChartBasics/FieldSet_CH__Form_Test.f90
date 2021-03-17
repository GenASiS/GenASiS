program FieldSet_CH__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( FieldSet_CH_Form ), allocatable :: &
    FSC
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    FSM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_CH__Form_Test', DimensionalityOption = '2D' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize_H &
         ( M, Periodic, iChart = 1 )
  call M % Show ( )
  call C % Show ( )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize ( M, 'Fields', iFieldSet = 1 ) 
  call FSC % Initialize ( FSM, C, 'Fields' ) 

  deallocate ( FSC )
  deallocate ( FSM )
  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_CH__Form_Test
