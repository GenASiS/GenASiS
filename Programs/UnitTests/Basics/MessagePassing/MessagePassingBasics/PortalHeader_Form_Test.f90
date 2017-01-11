program PortalHeader_Form_Test

  use Display
  use Communicator_Form
  use PortalHeader_Form

  implicit none
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( PortalHeaderForm ) :: &
    PH

  allocate ( C )
  
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )

  call PH % Initialize ( [ 0, 2 ], [ 1, 3 ], [ 10, 20 ], [ 100, 200 ] )

  call Show ( 'PortalHeaderForm members' )
  call Show ( PH % nSources, 'nSources' )
  call Show ( PH % nTargets, 'nTargets' )
  call Show ( PH % nChunksFrom, 'nChunksFrom' )
  call Show ( PH % nChunksTo, 'nChunksTo' )
  call Show ( PH % Source, 'Source' )
  call Show ( PH % Target, 'Target' )

  deallocate ( C )

end program PortalHeader_Form_Test
