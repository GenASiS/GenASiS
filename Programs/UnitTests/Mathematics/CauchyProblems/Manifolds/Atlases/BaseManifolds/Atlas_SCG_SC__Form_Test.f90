program Atlas_SCG_SC__Form_Test

  !-- Atlas_SingleChartGrid_SymmetricCurvilinear_Form_Test

  use Basics
  use BaseManifolds

  implicit none

  type ( Atlas_SCG_SC_Form ), allocatable :: &
    A

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Atlas_SCG_SC__Form_Test', DimensionalityOption = '2D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace' )

  call A % Show ( )
  
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program Atlas_SCG_SC__Form_Test
