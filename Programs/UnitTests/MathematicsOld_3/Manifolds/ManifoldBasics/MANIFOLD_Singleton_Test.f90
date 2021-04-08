program MANIFOLD_Singleton_Test

  use Basics
  use ManifoldBasics

  implicit none

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize ( 'MANIFOLD_Singleton_Test' )

  call Show ( MANIFOLD % MAX_DIMENSIONS, 'MANIFOLD % MAX_DIMENSIONS', &
              nLeadingLinesOption = 2 )
  call Show ( MANIFOLD % MAX_CHARTS,     'MANIFOLD % MAX_CHARTS' )
  call Show ( MANIFOLD % MAX_FIELD_SETS, 'MANIFOLD % MAX_FIELD_SETS' )
  call Show ( MANIFOLD % MAX_STREAMS,    'MANIFOLD % MAX_STREAMS', &
              nTrailingLinesOption = 2 )
  
  deallocate ( PROGRAM_HEADER )

end program MANIFOLD_Singleton_Test
