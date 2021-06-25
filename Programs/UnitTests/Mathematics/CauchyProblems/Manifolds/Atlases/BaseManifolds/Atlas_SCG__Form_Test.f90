program Atlas_SCG__Form_Test

  !-- Atlas_SingleChartGrid__Form_Test

  use Basics
  use Charts
  use BaseManifolds

  implicit none

  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  type ( Atlas_SCG_Form ), allocatable :: &
    Base, &
    Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Atlas_SCG__Form_Test', DimensionalityOption = '2D_1D' )

  allocate ( Base )
  call Base % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'Base', &
           iDimensionalityOption = 1 )

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( Fiber )
  call Fiber % Initialize &
         ( SpacingOption = [ 'GEOMETRIC' ], &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           NameOption = 'Fiber', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           MinCoordinateOption = [ MinEnergy ], &
           MaxCoordinateOption = [ MaxEnergy ], &
           ScaleOption = [ MinWidthEnergy ], &
           nCellsOption = [ 16 ], &
           nGhostLayersOption = [ 0 ], &
           iDimensionalityOption = 2 )

  call Base % Show ( )
  call Fiber % Show ( )

  deallocate ( Fiber )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program Atlas_SCG__Form_Test
