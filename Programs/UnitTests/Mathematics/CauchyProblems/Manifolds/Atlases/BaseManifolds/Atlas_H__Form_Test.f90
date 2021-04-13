program Atlas_H__Form_Test

  !-- Atlas_Header__Form_Test

  use Basics
  use Charts
  use BaseManifolds

  implicit none

  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  type ( Atlas_H_Form ), allocatable :: &
    A_Base, &
    A_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Atlas_H__Form_Test', DimensionalityOption = '2D_1D' )

  allocate ( A_Base )
  call A_Base % Initialize_H ( NameOption = 'Base' )

  allocate ( Grid_S_Form :: A_Base % Chart ( 1 ) % Element )
  select type ( G_Base  =>  A_Base % Chart ( 1 ) % Element )
  class is ( Grid_S_Form )
  call G_Base % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'G_Base', &
           PeriodicOption = [ .true., .true., .true. ], &
           iDimensionalityOption = 1 )
  end select !-- G_Base

  allocate ( A_Fiber )
  call A_Fiber % Initialize_H ( NameOption = 'Fiber' )

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( Grid_S_Form :: A_Fiber % Chart ( 1 ) % Element )
  select type ( G_Fiber  =>  A_Fiber % Chart ( 1 ) % Element )
  class is ( Grid_S_Form )
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
  end select !-- G_Fiber

  call A_Base % Show ( )
  call A_Fiber % Show ( )

  deallocate ( A_Fiber )
  deallocate ( A_Base )
  deallocate ( PROGRAM_HEADER )

end program Atlas_H__Form_Test
