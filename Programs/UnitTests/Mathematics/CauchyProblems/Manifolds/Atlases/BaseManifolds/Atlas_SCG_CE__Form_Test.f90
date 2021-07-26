program Atlas_SCG_CE__Form_Test

  !-- Atlas_SingleChartGrid_CentralExcision_Form_Test

  use Basics
  use BaseManifolds

  implicit none

  type ( Atlas_SCG_CE_Form ), allocatable :: &
    A, &
    A_SA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Atlas_SCG_CE__Form_Test', DimensionalityOption = '2D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusExcision = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace' )

  allocate ( A_SA )
  call A_SA % Initialize ( A )

  call A   % Show ( )
  call A_SA % Show ( )

  deallocate ( A_SA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program Atlas_SCG_CE__Form_Test
