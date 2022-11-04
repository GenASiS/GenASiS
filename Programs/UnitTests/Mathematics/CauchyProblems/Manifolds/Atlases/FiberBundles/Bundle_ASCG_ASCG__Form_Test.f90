program Bundle_ASCG_ASCG__Form_Test

  !-- Atlas_SingleChartGrid__Form_Test

  use Basics
  use Charts
  use BaseManifolds
  use FiberBundles

  implicit none

  
  integer ( KDI ) :: &
    iPF, &  !-- iProcessFiber
    iR, &   !-- iRank
    nProcesses, &
    nProcessesBase, &
    nProcessesFiber
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
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  Communicator     =>  PROGRAM_HEADER % Communicator
  nProcesses       =   Communicator % Size
  nProcessesFiber  =   1

  call PROGRAM_HEADER % GetParameter ( nProcessesFiber, 'nProcessesFiber' )
  nProcessesBase  =  nProcesses / nProcessesFiber

  if ( nProcessesBase * nProcessesFiber  /=  nProcesses ) then
    call Show ( 'nProcessesBase * nProcessesFiber must equal nProcesses', &
                CONSOLE % ERROR )
    call Show ( nProcesses,      'nProcesses',      CONSOLE % ERROR )
    call Show ( nProcessesBase,  'nProcessesBase',  CONSOLE % ERROR )
    call Show ( nProcessesFiber, 'nProcessesFiber', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )
  end if

  call Show ( nProcesses,      'nProcesses'      )
  call Show ( nProcessesBase,  'nProcessesBase'  )
  call Show ( nProcessesFiber, 'nProcessesFiber' )

  !-- nProcessesFiber copies of Base

  do iPF  =  1, nProcessesFiber
    allocate ( Rank, &
!               source =  [ ( iR, iR = iPF - 1, nProcesses - 1, &
!                                      nProcessesFiber ) ] )
               source =  [ ( iR, iR = ( iPF - 1 ) * nProcessesBase, &
                                      iPF * nProcessesBase  -  1 ) ] )
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
  end do !--iPF

  ! !-- Bundle

  !      MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  !      MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  ! MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  ! allocate ( Bundle )
  ! call Bundle % Initialize &
  !        ( Base, &
  !          SpacingOption = [ 'GEOMETRIC' ], &
  !          CoordinateLabelOption = [ 'E' ], &
  !          CoordinateSystemOption = 'SPHERICAL', &
  !          NameOption = 'Fiber', &
  !          CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
  !          MinCoordinateOption = [ MinEnergy ], &
  !          MaxCoordinateOption = [ MaxEnergy ], &
  !          ScaleOption = [ MinWidthEnergy ], &
  !          nCellsOption = [ 16 ], &
  !          nGhostLayersOption = [ 0 ] )

  call Base % Show ( )
  ! call Bundle % Show ( )

  ! deallocate ( Bundle )
  deallocate ( Base )
  deallocate ( CommunicatorBase )
  deallocate ( PROGRAM_HEADER )

end program Bundle_ASCG_ASCG__Form_Test
