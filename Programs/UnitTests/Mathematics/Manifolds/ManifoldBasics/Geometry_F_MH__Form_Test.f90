program Geometry_F_MH__Form_Test

  use Basics
  use ManifoldBasics

  implicit none

  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( Geometry_F_MH_Form ), allocatable :: &
    GM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F_MH__Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( GM )
  call GM % Initialize_F ( M ) 

  call M % Show ( )
  call Show ( M % nCharts,    'nCharts',    M % IGNORABILITY )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call GM % Show ( )
  call Show ( GM % nStreams,  'nStreams',  GM % IGNORABILITY )

  deallocate ( GM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_MH__Form_Test
