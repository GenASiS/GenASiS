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

  allocate ( A_Base % Chart ( 1 ) % Element )
  associate ( C_Base  =>  A_Base % Chart ( 1 ) % Element )
  call C_Base % Initialize_H &
         ( NameOption = 'C_Base', &
           iDimensionalityOption = 1 )
  end associate !-- C_Base

  allocate ( A_Fiber )
  call A_Fiber % Initialize_H ( NameOption = 'Fiber' )

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( A_Fiber % Chart ( 1 ) % Element )
  associate ( C_Fiber  =>  A_Fiber % Chart ( 1 ) % Element )
  call C_Fiber % Initialize_H &
         ( CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           NameOption = 'C_Fiber', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           iDimensionalityOption = 2 )
  end associate !-- C_Fiber

  call A_Base % Show ( )
  call A_Fiber % Show ( )

  deallocate ( A_Fiber )
  deallocate ( A_Base )
  deallocate ( PROGRAM_HEADER )

end program Atlas_H__Form_Test
