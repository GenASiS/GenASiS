program Chart_BH__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics
  use BaseCharts

  implicit none

  real ( KDR ) :: &
    MinEnergy, &
    MaxEnergy, &
    MinWidthEnergy
  logical ( KDL ), dimension ( 3 ) :: &
    IsPeriodic
  type ( Manifold_H_Form ), allocatable :: &
    Base, &
    Fiber
  type ( FieldSet_MH_Form ), allocatable :: &
    GM_Base, &
    GM_Fiber
  type ( FieldSet_CH_Form ), allocatable :: &
    GC_Base, &
    GC_Fiber
  type ( Geometry_F_Form ), allocatable :: &
    G_Base, &
    G_Fiber
  type ( Chart_BH_Form ), allocatable :: &
    C_Base, &
    C_Fiber

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Chart_BH__Form_Test', DimensionalityOption = '2D_1D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  !-- Base

  IsPeriodic  =  .true.

  allocate ( Base )
  allocate ( GM_Base )
  allocate ( C_Base )
  allocate ( GC_Base )
  allocate ( G_Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call GM_Base % Initialize &
         ( Base, 'Geometry' ) 
  call C_Base % Initialize_BH &
         ( Base, IsPeriodic, iChart = 1 )
  call GC_Base % Initialize &
         ( GM_Base, C_Base, 'Geometry' ) 
  call G_Base % Initialize &
         ( GC_Base, nValues = C_Base % nValues, NameOption = 'Geometry' )
  call C_Base % ComputeGeometry ( G_Base )

  call Base % Show ( )
  call C_Base % Show ( )

  !-- Fiber

  IsPeriodic  =  .false.

       MinEnergy  =    0.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
       MaxEnergy  =  100.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  MinWidthEnergy  =    0.1_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  allocate ( Fiber )
  allocate ( GM_Fiber )
  allocate ( C_Fiber )
  allocate ( GC_Fiber )
  allocate ( G_Fiber )
  call Fiber % Initialize &
         ( 'Fiber', iDimensionalityOption = 2 )
  call GM_Fiber % Initialize &
         ( Fiber, 'Geometry' ) 
  call C_Fiber % Initialize_BH &
         ( Fiber, IsPeriodic, iChart = 1, &
           SpacingOption = [ 'GEOMETRIC' ], &
           CoordinateLabelOption = [ 'E' ], &
           CoordinateSystemOption = 'SPHERICAL', &
           CoordinateUnitOption = [ UNIT % MEGA_ELECTRON_VOLT ], &
           MinCoordinateOption = [ MinEnergy ], &
           MaxCoordinateOption = [ MaxEnergy ], &
           ScaleOption = [ MinWidthEnergy ], &
           nCellsOption = [ 16 ], &
           nGhostLayersOption = [ 0 ] )
  call GC_Fiber % Initialize &
         ( GM_Fiber, C_Fiber, 'Geometry' ) 
  call G_Fiber % Initialize &
         ( GC_Fiber, nValues = C_Fiber % nValues, NameOption = 'Geometry' )
  call C_Fiber % ComputeGeometry ( G_Fiber )

  call Fiber % Show ( )
  call C_Fiber % Show ( )

  deallocate ( G_Fiber )
  deallocate ( GC_Fiber )
  deallocate ( C_Fiber )
  deallocate ( GM_Fiber )
  deallocate ( Fiber )
  deallocate ( G_Base )
  deallocate ( GC_Base )
  deallocate ( C_Base )
  deallocate ( GM_Base )
  deallocate ( Base )
  deallocate ( PROGRAM_HEADER )

end program Chart_BH__Form_Test
