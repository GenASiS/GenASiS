program FieldSet_CH__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  integer ( KDI ) :: &
    nFields = 5
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
         ( M, 'Global', Periodic )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize ( M, 'Fields', nFields ) 
  call FSC % Initialize ( C, FSM ) 

  call M % Show ( )
  call Show ( M % nCharts,    'nCharts',    M % IGNORABILITY )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call FSM % Show ( )
  call Show ( FSM % nStreams,  'nStreams',  FSM % IGNORABILITY )

  call C % Show ( )
  call Show ( C % nFieldSets, 'nFieldSets', C % IGNORABILITY )
  call Show ( C % nStreams,   'nStreams',   C % IGNORABILITY )

  call FSC % Show ( )
  call Show ( FSC % nStreams,  'nStreams',  FSC % IGNORABILITY )

  deallocate ( FSC )
  deallocate ( FSM )
  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_CH__Form_Test
