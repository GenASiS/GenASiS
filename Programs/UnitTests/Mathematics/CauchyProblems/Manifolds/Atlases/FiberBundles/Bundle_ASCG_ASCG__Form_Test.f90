program Bundle_ASCG_ASCG__Form_Test

  !-- Atlas_SingleChartGrid__Form_Test

  use Basics
  use Charts
  use BaseManifolds
  use FiberBundles

  implicit none

  
  integer ( KDI ) :: &
    iCB, &  !-- iCopyBase
    iR, &   !-- iRank
    nProcesses, &
    nProcessesBase, &
    nCopiesBase, &  !-- for distribution of operations on the base manifold
                    !   for different fiber bins
    nFiberBins
  integer ( KDI ), dimension ( : ), allocatable :: &
    Rank
  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  type ( CommunicatorForm ), allocatable :: &
    CommunicatorBase
  type ( CommunicatorForm ), pointer :: &
    Communicator
  type ( Atlas_SCG_Form ), allocatable :: &
    Base
  type ( Bundle_ASCG_ASCG_Form ), allocatable :: &
    Bundle

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Bundle_ASCG_ASCG__Form_Test', DimensionalityOption = '2D_1D' )
  call CONSOLE % SetVerbosity ( 'INFO_3' )

  Communicator  =>  PROGRAM_HEADER % Communicator
  nProcesses    =   Communicator % Size
  nCopiesBase   =   1

  call PROGRAM_HEADER % GetParameter ( nCopiesBase, 'nCopiesBase' )
  nProcessesBase  =  nProcesses / nCopiesBase

  call Show ( nProcesses,     'nProcesses'      )
  call Show ( nProcessesBase, 'nProcessesBase'  )
  call Show ( nCopiesBase,    'nCopiesBase' )

  do iCB  =  1, nCopiesBase
    allocate ( Rank, &
!               source =  [ ( iR, iR = iCB - 1, nProcesses - 1, &
!                                      nCopiesBase ) ] )
               source =  [ ( iR, iR = ( iCB - 1 ) * nProcessesBase, &
                                      iCB * nProcessesBase  -  1 ) ] )
    if ( any ( Communicator % Rank  ==  Rank ) ) then
      allocate ( CommunicatorBase )
      allocate ( Base )
      call CommunicatorBase % Initialize &
             ( Communicator, Rank, NameOption = 'CommunicatorBase' )
      call Base % Initialize &
             ( CommunicatorOption = CommunicatorBase, &
               NameOption = 'Base', &
               iDimensionalityOption = 1 )
    end if
    deallocate ( Rank )
  end do !--iCB

  !-- Bundle

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  nFiberBins  =  16
  call PROGRAM_HEADER % GetParameter ( nFiberBins, 'nFiberBins' )

  allocate ( Bundle )
  call Bundle % Initialize &
         ( Base, &
           SpacingOption = [ 'GEOMETRIC' ], &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           NameOption = 'Fiber', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           MinCoordinateOption = [ MinEnergy ], &
           MaxCoordinateOption = [ MaxEnergy ], &
           ScaleOption = [ MinWidthEnergy ], &
           nCellsOption = [ nFiberBins ], &
           nGhostLayersOption = [ 0 ] )

  call Base % Show ( )
  call Bundle % Show ( )

  deallocate ( Bundle )
  deallocate ( Base )
  deallocate ( CommunicatorBase )
  deallocate ( PROGRAM_HEADER )

end program Bundle_ASCG_ASCG__Form_Test
