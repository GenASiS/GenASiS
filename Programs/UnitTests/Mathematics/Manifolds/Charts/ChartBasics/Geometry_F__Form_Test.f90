program Geometry_F__Form_Test

  !-- Geometry_Flat_Form

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    GM
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( FieldSet_CH_Form ), allocatable :: &
    GC
  type ( Geometry_F_Form ), allocatable :: &
    G

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'Geometry_F__Form_Test' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize_H &
         ( M, 'Global', Periodic )

  allocate ( GM )
  allocate ( GC )
  call GM % Initialize ( M, 'Geometry' ) 
  call GC % Initialize ( C, GM, nFields = 19 ) 

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  allocate ( G )
  call G % Initialize ( GC, nValues = 10 )
  deallocate ( G )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call M % Show ( )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call GM % Show ( )
  call Show ( GM % nStreams,  'nStreams',  GM % IGNORABILITY )

  call C % Show ( )
  call GC % Show ( )

  deallocate ( GC )
  deallocate ( GM )
  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F__Form_Test
