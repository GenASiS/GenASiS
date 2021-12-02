program Coarsening_C__Form_Test

  !-- Coarsening_C__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( Coarsening_C_Form ), allocatable :: &
    C

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Coarsening_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace', &
           nCellsPolarOption = 32 )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( C )
  call C % Initialize ( G )
  call S % AddFieldSet ( C )

  call A % Show ( )
  call G % Show ( )
  call C % Show ( )
  call S % Show ( )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  deallocate ( C )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Coarsening_C__Form_Test
