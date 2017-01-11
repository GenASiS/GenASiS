program ATLAS_Singleton_Test

  use Basics
  use ATLAS_Singleton

  implicit none

  character ( LDF ) :: &
    Name = 'ATLAS_Singleton_Test'

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize ( Name )

  call Show ( ATLAS % MAX_DIMENSIONS, 'ATLAS % MAX_DIMENSIONS', &
              nLeadingLinesOption = 2 )
  call Show ( ATLAS % MAX_CHARTS, 'ATLAS % MAX_CHARTS' )
  call Show ( ATLAS % MAX_CHARTS, 'ATLAS % MAX_FIELDS' )
  call Show ( ATLAS % MAX_STREAMS, 'ATLAS % MAX_STREAMS', &
              nTrailingLinesOption = 2 )
  
  deallocate ( PROGRAM_HEADER )

end program ATLAS_Singleton_Test
