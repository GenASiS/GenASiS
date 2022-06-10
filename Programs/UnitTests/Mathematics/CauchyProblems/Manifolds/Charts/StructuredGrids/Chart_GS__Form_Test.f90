program Chart_GS__Form_Test

  !-- Chart_GridStructured__Form_Test

  use Basics
  use StructuredGrids

  implicit none

  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  type ( Chart_GS_Form ), allocatable :: &
    C_Base, &
    C_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Chart_GS__Form_Test', DimensionalityOption = '2D_1D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( C_Base )
  call C_Base % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'C_Base', &
           iDimensionalityOption = 1 )

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( C_Fiber )
  call C_Fiber % Initialize &
         ( SpacingOption = [ 'GEOMETRIC' ], &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           NameOption = 'C_Fiber', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           MinCoordinateOption = [ MinEnergy ], &
           MaxCoordinateOption = [ MaxEnergy ], &
           ScaleOption = [ MinWidthEnergy ], &
           nCellsOption = [ 16 ], &
           nGhostLayersOption = [ 0 ], &
           iDimensionalityOption = 2 )

  call C_Base % Show ( )
  call C_Fiber % Show ( )

  deallocate ( C_Fiber )
  deallocate ( C_Base )
  deallocate ( PROGRAM_HEADER )

end program Chart_GS__Form_Test
