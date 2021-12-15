program Coarsening_C__Form_Test

  !-- Coarsening_C__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none

  real ( KDR ), dimension ( :, :, : ), pointer :: &
    FV_3D, &
    Th_3D, &
    Ph_3D
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( FieldSetForm ), allocatable :: &
    F
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

  allocate ( F )
  call F % Initialize &
         ( A, &
           FieldOption = [ 'Field' ], &
           NameOption = 'Field', &
           nFieldsOption = 1 )
  call S % AddFieldSet ( F )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( C )
  call C % Initialize ( G )
  call S % AddFieldSet ( C )

  call A % Show ( )
  call G % Show ( )

  call C % Show ( )
  call Show ( C % nBlocksCoarsen, 'nBlocksCoarsen' )
  call Show ( C % iRadius, 'iRadius' )
  call Show ( C % iTheta, 'iTheta' )
  call Show ( C % iPhi, 'iPhi' )

  call S % Show ( )

  call A % Chart_GS % SetFieldPointer &
         ( F % Storage_GS % Value ( :, 1 ), FV_3D )
  call A % Chart_GS % SetFieldPointer &
         ( G % Storage_GS % Value ( :, G % CENTER_U_2 ), Th_3D )
  call A % Chart_GS % SetFieldPointer &
         ( G % Storage_GS % Value ( :, G % CENTER_U_3 ), Ph_3D )
  select case ( A % Chart_GS % nDimensions )
  case ( 2 )
    FV_3D  =  sin ( Th_3D )
  case ( 3 )
    FV_3D  =  sin ( Th_3D )  *  sin ( Ph_3D )
  end select !-- nDimensions

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  call C % Compute ( F )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  deallocate ( C )
  deallocate ( G )
  deallocate ( F )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Coarsening_C__Form_Test
