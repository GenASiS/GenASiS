program Grid_S__Form_Test

  !-- Grid_Structured__Form_Test

  use Basics
  use StructuredGrids

  implicit none

  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  type ( Grid_S_Form ), allocatable :: &
    G_Base, &
    G_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Grid_S__Form_Test', DimensionalityOption = '2D_1D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( G_Base )
  call G_Base % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'G_Base', &
           PeriodicOption = [ .true., .true., .true. ], &
           iDimensionalityOption = 1 )

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( G_Fiber )
  call G_Fiber % Initialize &
         ( SpacingOption = [ 'GEOMETRIC' ], &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           NameOption = 'G_Fiber', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           MinCoordinateOption = [ MinEnergy ], &
           MaxCoordinateOption = [ MaxEnergy ], &
           ScaleOption = [ MinWidthEnergy ], &
           nCellsOption = [ 16 ], &
           nGhostLayersOption = [ 0 ], &
           iDimensionalityOption = 2 )

  call G_Base % Show ( )
  call G_Fiber % Show ( )

  deallocate ( G_Fiber )
  deallocate ( G_Base )
  deallocate ( PROGRAM_HEADER )

end program Grid_S__Form_Test
