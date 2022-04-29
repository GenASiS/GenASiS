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
  type ( Atlas_SCG_CE_Form ), allocatable :: &
    A_CE
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A_CC
  type ( FieldSetForm ), allocatable :: &
    F_CE, &
    F_CC
  type ( StreamForm ), allocatable :: &
    S_CE, &
    S_CC
  type ( Geometry_F_Form ), allocatable :: &
    G_CE, &
    G_CC
  type ( Coarsening_C_Form ), allocatable :: &
    C_CE, &
    C_CC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Coarsening_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A_CE, A_CC )
  call A_CE % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusExcision = 0.45_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace_CE', &
           nCellsPolarOption = 32 )
  call A_CC % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace_CC', &
           nCellsPolarOption = 32 )

  allocate ( S_CE, S_CC )
  call S_CE % Initialize ( A_CE, GIS )
  call S_CC % Initialize ( A_CC, GIS )

  allocate ( F_CE, F_CC )
  call F_CE % Initialize &
         ( A_CE, &
           FieldOption = [ 'Field_CE' ], &
           NameOption = 'Field_CE', &
           nFieldsOption = 1 )
  call F_CC % Initialize &
         ( A_CC, &
           FieldOption = [ 'Field_CC' ], &
           NameOption = 'Field_CC', &
           nFieldsOption = 1 )
  call S_CE % AddFieldSet ( F_CE )
  call S_CC % AddFieldSet ( F_CC )

  allocate ( G_CE, G_CC )
  call G_CE % Initialize ( A_CE )
  call G_CC % Initialize ( A_CC )
  call G_CE % SetStream ( S_CE )
  call G_CC % SetStream ( S_CC )

  allocate ( C_CE, C_CC )
  call C_CE % Initialize ( G_CE, NameOption = 'Coarsening_CE' )
  call C_CC % Initialize ( G_CC, NameOption = 'Coarsening_CC' )
  call S_CE % AddFieldSet ( C_CE )
  call S_CC % AddFieldSet ( C_CC )

  call A_CE % Show ( )
  call A_CC % Show ( )

  call G_CE % Show ( )
  call G_CC % Show ( )

  call C_CE % Show ( )
  call C_CC % Show ( )
  call Show ( C_CE % nBlocksCoarsen, 'CE nBlocksCoarsen' )
  call Show ( C_CE % iRadius, 'CE iRadius' )
  call Show ( C_CE % iTheta, 'CE iTheta' )
  call Show ( C_CE % iPhi, 'CE iPhi' )
  call Show ( C_CC % nBlocksCoarsen, 'CC nBlocksCoarsen' )
  call Show ( C_CC % iRadius, 'CC iRadius' )
  call Show ( C_CC % iTheta, 'CC iTheta' )
  call Show ( C_CC % iPhi, 'CC iPhi' )

  call S_CE % Show ( )
  call S_CC % Show ( )

  call A_CC % Chart_GS % SetFieldPointer &
         ( F_CC % Storage_GS % Value ( :, 1 ), FV_3D )
  call A_CC % Chart_GS % SetFieldPointer &
         ( G_CC % Storage_GS % Value ( :, G_CC % CENTER_U_2 ), Th_3D )
  call A_CC % Chart_GS % SetFieldPointer &
         ( G_CC % Storage_GS % Value ( :, G_CC % CENTER_U_3 ), Ph_3D )
  select case ( A_CC % Chart_GS % nDimensions )
  case ( 2 )
    FV_3D  =  sin ( Th_3D )
  case ( 3 )
    FV_3D  =  sin ( Th_3D )  *  sin ( Ph_3D )
  end select !-- nDimensions

  call A_CE % Chart_GS % SetFieldPointer &
         ( F_CE % Storage_GS % Value ( :, 1 ), FV_3D )
  call A_CE % Chart_GS % SetFieldPointer &
         ( G_CE % Storage_GS % Value ( :, G_CE % CENTER_U_2 ), Th_3D )
  call A_CE % Chart_GS % SetFieldPointer &
         ( G_CE % Storage_GS % Value ( :, G_CE % CENTER_U_3 ), Ph_3D )
  select case ( A_CE % Chart_GS % nDimensions )
  case ( 2 )
    FV_3D  =  sin ( Th_3D )
  case ( 3 )
    FV_3D  =  sin ( Th_3D )  *  sin ( Ph_3D )
  end select !-- nDimensions

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S_CE % Write ( )
  call S_CC % Write ( )
  call GIS % Close ( )

  call C_CE % Compute ( F_CE )
  call C_CC % Compute ( F_CC )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S_CE % Write ( )
  call S_CC % Write ( )
  call GIS % Close ( )

  deallocate ( C_CC, C_CE )
  deallocate ( G_CC, G_CE )
  deallocate ( F_CC, F_CE )
  deallocate ( S_CC, S_CE )
  deallocate ( A_CC, A_CE )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Coarsening_C__Form_Test
