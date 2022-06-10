program Atlas_SCG_CC__Form_Test

  !-- Atlas_SingleChartGrid_CentralCore_Form_Test

  use Basics
  use BaseManifolds

  implicit none

  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A, &
    A_SA, &  !-- SphericalAverage
    A_AA     !-- AzimuthalAverage

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Atlas_SCG_CC__Form_Test', DimensionalityOption = '2D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace' )

  if ( A % Chart_GS_CC % nDimensions  >  1 ) then
    allocate ( A_SA )
    call A_SA % Initialize ( A, nDimensions = 1 )
  end if

  if ( A % Chart_GS_CC % nDimensions  >  2 ) then
    allocate ( A_AA )
    call A_AA % Initialize ( A, nDimensions = 2 )
  end if

  call A % Show ( )
  if ( allocated ( A_SA ) ) &
    call A_SA % Show ( )
  if ( allocated ( A_AA ) ) &
    call A_AA % Show ( )
  
  if ( allocated ( A_AA ) ) &
    deallocate ( A_AA )
  if ( allocated ( A_SA ) ) &
    deallocate ( A_SA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program Atlas_SCG_CC__Form_Test
