program CurrentSet_A__Form_Test

  !-- CurrentSet_Atlas__Form_Test

  use Basics
  use Manifolds
  use Streams
  use CurrentSets

  implicit none

  type ( MeasuredValueForm ) :: &
    DensityUnit
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    Velocity_U_Unit
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( CurrentSet_A_Form ), allocatable :: &
    CSA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_F_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  Velocity_U_Unit  =  UNIT % SPEED_MKS
      DensityUnit  =  UNIT % MASS_DENSITY_MKS

  allocate ( CSA )
  call CSA % Initialize &
         ( A, &
           Velocity_U_Unit, &
           DensityUnitOption = DensityUnit )
  call CSA % SetStream ( SA )

  call   A % Show ( )
  call CSA % Show ( )
  call  SA % Show ( )

  deallocate ( CSA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program CurrentSet_A__Form_Test
