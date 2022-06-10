program Fluid_P_I__Form_Test

  !-- Fluid_Perfect_Ideal__Form_Test

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( Gravitation_G_Form ), allocatable :: &
    G
  type ( Units_F_Form ), dimension ( : ), allocatable :: &
    U
  type ( Fluid_P_I_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Fluid_P_I__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( U ( 1 ) )
  call U ( 1 ) % Initialize ( TypeOption = 'MKS' )

  allocate ( F )
  call F % Initialize ( G, U )
  call F % SetStream ( S )

  call A       % Show ( )
  call U ( 1 ) % Show ( )
  call F       % Show ( )
  call S       % Show ( )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  deallocate ( F )
  deallocate ( U )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Fluid_P_I__Form_Test
