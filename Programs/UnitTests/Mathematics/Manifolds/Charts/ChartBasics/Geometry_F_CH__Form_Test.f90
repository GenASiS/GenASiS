program Geometry_F_CH__Form_Test

  !-- Geometry_Flat_ChartHeader_Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( Geometry_F_MH_Form ), allocatable :: &
    GM
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( Geometry_F_CH_Form ), allocatable :: &
    GC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'Geometry_F_CH__Form_Test' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( GM )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call GM % Initialize_F ( M ) 

  allocate ( C )
  allocate ( GC )
  call C % Initialize_H &
        ( M, 'Global', Periodic )
  call GC % Initialize ( C, GM ) 

  call M % Show ( )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call GM % Show ( )
  call Show ( GM % nStreams,  'nStreams',  GM % IGNORABILITY )

  call C % Show ( )
  call GC % Show ( )

  deallocate ( GC )
  deallocate ( C )
  deallocate ( GM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_CH__Form_Test
