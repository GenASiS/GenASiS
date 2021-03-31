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
  real ( KDR ), dimension ( :, :, : ), pointer :: &
    Edge_I_3D, &
    Width_3D, &
    Center_3D, &
    Area_I_3D, &
    Volume_3D, &
    Metric_F_DD_11_3D, &
    Metric_F_DD_22_3D, &
    Metric_F_DD_33_3D, &
    Metric_F_UU_11_3D, &
    Metric_F_UU_22_3D, &
    Metric_F_UU_33_3D
  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( Manifold_H_Form ), allocatable :: &
    Base, &
    Fiber
  type ( Geometry_F_MH_Form ), allocatable :: &
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

  Periodic  =  .true.

  allocate ( Base )
  allocate ( GM_Base )
  allocate ( C_Base )
  allocate ( GC_Base )
  allocate ( G_Base )
  call Base % Initialize &
         ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
           iDimensionalityOption = 1 )
  call GM_Base % Initialize &
         ( Base, NameOption = 'GeometryBase' ) 
  call C_Base % Initialize &
         ( Base, 'Global', Periodic )
  call GC_Base % Initialize &
         ( C_Base, GM_Base ) 
  call G_Base % Initialize &
         ( GC_Base, nValues = C_Base % nValues )
  call C_Base % ComputeGeometry ( G_Base )

  !-- Fiber

  Periodic  =  .false.

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
         ( Fiber, NameOption = 'GeometryFiber' ) 
  call C_Fiber % Initialize &
         ( Fiber, 'Global', Periodic, &
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
         ( C_Fiber, GM_Fiber ) 
  call G_Fiber % Initialize &
         ( GC_Fiber, nValues = C_Fiber % nValues )
  call C_Fiber % ComputeGeometry ( G_Fiber )

  !-- Display and cleanup

  call Base % Show ( )
  call Show ( Base % nCharts,    'nCharts',    Base % IGNORABILITY )
  call Show ( Base % nFieldSets, 'nFieldSets', Base % IGNORABILITY )
  call Show ( Base % nStreams,   'nStreams',   Base % IGNORABILITY )

  call C_Base % Show ( )

  call ShowProper ( C_Base )
  call ShowGeometry ( G_Base )

  call Fiber % Show ( )
  call Show ( Fiber % nCharts,    'nCharts',    Fiber % IGNORABILITY )
  call Show ( Fiber % nFieldSets, 'nFieldSets', Fiber % IGNORABILITY )
  call Show ( Fiber % nStreams,   'nStreams',   Fiber % IGNORABILITY )

  call C_Fiber % Show ( )

  call ShowProper ( C_Fiber )
  call ShowGeometry ( G_Fiber )

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

contains


  subroutine ShowProper ( C )

    class ( Chart_BH_Form ), intent ( inout ) :: &
      C

    logical ( KDL ), dimension ( :, :, : ), pointer :: &
      PC

    call Show ( 'IsProper' )

    associate &
      ( iaF  =>  C % iaFirst, &
        iaL  =>  C % iaLast )
    PC ( iaF ( 1 ) : iaL ( 1 ), &
         iaF ( 2 ) : iaL ( 2 ), &
         iaF ( 3 ) : iaL ( 3 ) )  &
      =>  C % ProperCell
    end associate !-- iaF

    call Show ( 1, 'iDimension' )
    call Show ( PC ( :, 1, 1 ), 'ProperCell' )

    if ( C % nDimensions > 1 ) then
      call Show ( 2, 'iDimension' )
      call Show ( PC ( 1, :, 1 ), 'ProperCell' )
    end if

    if ( C % nDimensions > 2 ) then
      call Show ( 3, 'iDimension' )
      call Show ( PC ( 1, 1, : ), 'ProperCell' )
    end if

    nullify ( PC )

  end subroutine ShowProper


  subroutine ShowGeometry ( G )

    type ( Geometry_F_Form ), intent ( in ) :: &
      G

    call Show ( G % Name )

    select type ( C  =>  G % Geometry_C % Chart )
    class is ( Chart_BH_Form )

    call Show ( 1, 'iDimension' )
    call SetGeometryPointers ( G, iD = 1 )
    call Show ( Edge_I_3D ( :, 1, 1 ), G % Unit ( G % EDGE_I_U ( 1 ) ), &
                'Edge_I' )
    call Show (  Width_3D ( :, 1, 1 ), G % Unit ( G % WIDTH_U  ( 1 ) ), &
                 'Width' )
    call Show ( Center_3D ( :, 1, 1 ), G % Unit ( G % CENTER_U ( 1 ) ), &
                'Center' )
    call Show ( Area_I_3D ( :, 1, 1 ), G % Unit ( G % AREA_I_D ( 1 ) ), &
                'Area_I' )
    call Show ( Volume_3D ( :, 1, 1 ), G % Unit ( G % VOLUME ), &
                'Volume' )
    call Show ( Metric_F_DD_11_3D ( :, 1, 1 ), &
                G % Unit ( G % METRIC_F_DD_11 ), 'Metric_F_DD_11' )
    call Show ( Metric_F_DD_22_3D ( :, 1, 1 ), &
                G % Unit ( G % METRIC_F_DD_22 ), 'Metric_F_DD_22' )
    call Show ( Metric_F_DD_33_3D ( :, 1, 1 ), &
                G % Unit ( G % METRIC_F_DD_33 ), 'Metric_F_DD_33' )
    call Show ( Metric_F_UU_11_3D ( :, 1, 1 ), &
                G % Unit ( G % METRIC_F_UU_11 ), 'Metric_F_UU_11' )
    call Show ( Metric_F_UU_22_3D ( :, 1, 1 ), &
                G % Unit ( G % METRIC_F_UU_22 ), 'Metric_F_UU_22' )
    call Show ( Metric_F_UU_33_3D ( :, 1, 1 ), &
                G % Unit ( G % METRIC_F_UU_33 ), 'Metric_F_UU_33' )

    if ( C % nDimensions > 1 ) then
      call Show ( 2, 'iDimension' )
      call SetGeometryPointers ( G, iD = 2 )
      call Show ( Edge_I_3D ( 1, :, 1 ), G % Unit ( G % EDGE_I_U ( 2 ) ), &
                  'Edge_I' )
      call Show (  Width_3D ( 1, :, 1 ), G % Unit ( G % WIDTH_U  ( 2 ) ), &
                   'Width' )
      call Show ( Center_3D ( 1, :, 1 ), G % Unit ( G % CENTER_U ( 2 ) ), &
                  'Center' )
      call Show ( Area_I_3D ( 1, :, 1 ), G % Unit ( G % AREA_I_D ( 2 ) ), &
                  'Area_I' )
      call Show ( Volume_3D ( 1, :, 1 ), G % Unit ( G % VOLUME ), &
                  'Volume' )
      call Show ( Metric_F_DD_11_3D ( 1, :, 1 ), &
                  G % Unit ( G % METRIC_F_DD_11 ), 'Metric_F_DD_11' )
      call Show ( Metric_F_DD_22_3D ( 1, :, 1 ), &
                  G % Unit ( G % METRIC_F_DD_22 ), 'Metric_F_DD_22' )
      call Show ( Metric_F_DD_33_3D ( 1, :, 1 ), &
                  G % Unit ( G % METRIC_F_DD_33 ), 'Metric_F_DD_33' )
      call Show ( Metric_F_UU_11_3D ( 1, :, 1 ), &
                  G % Unit ( G % METRIC_F_UU_11 ), 'Metric_F_UU_11' )
      call Show ( Metric_F_UU_22_3D ( 1, :, 1 ), &
                  G % Unit ( G % METRIC_F_UU_22 ), 'Metric_F_UU_22' )
      call Show ( Metric_F_UU_33_3D ( 1, :, 1 ), &
                  G % Unit ( G % METRIC_F_UU_33 ), 'Metric_F_UU_33' )
    end if

    if ( C % nDimensions > 2 ) then
      call Show ( 3, 'iDimension' )
      call SetGeometryPointers ( G, iD = 3 )
      call Show ( Edge_I_3D ( 1, 1, : ), G % Unit ( G % EDGE_I_U ( 3 ) ), &
                  'Edge_I' )
      call Show (  Width_3D ( 1, 1, : ), G % Unit ( G % WIDTH_U  ( 3 ) ), &
                   'Width' )
      call Show ( Center_3D ( 1, 1, : ), G % Unit ( G % CENTER_U ( 3 ) ), &
                  'Center' )
      call Show ( Area_I_3D ( 1, 1, : ), G % Unit ( G % AREA_I_D ( 3 ) ), &
                  'Area_I' )
      call Show ( Volume_3D ( 1, 1, : ), G % Unit ( G % VOLUME ), &
                  'Volume' )
      call Show ( Metric_F_DD_11_3D ( 1, 1, : ), &
                  G % Unit ( G % METRIC_F_DD_11 ), 'Metric_F_DD_11' )
      call Show ( Metric_F_DD_22_3D ( 1, 1, : ), &
                  G % Unit ( G % METRIC_F_DD_22 ), 'Metric_F_DD_22' )
      call Show ( Metric_F_DD_33_3D ( 1, 1, : ), &
                  G % Unit ( G % METRIC_F_DD_33 ), 'Metric_F_DD_33' )
      call Show ( Metric_F_UU_11_3D ( 1, 1, : ), &
                  G % Unit ( G % METRIC_F_UU_11 ), 'Metric_F_UU_11' )
      call Show ( Metric_F_UU_22_3D ( 1, 1, : ), &
                  G % Unit ( G % METRIC_F_UU_22 ), 'Metric_F_UU_22' )
      call Show ( Metric_F_UU_33_3D ( 1, 1, : ), &
                  G % Unit ( G % METRIC_F_UU_33 ), 'Metric_F_UU_33' )
    end if


    end select !-- C

  end subroutine ShowGeometry


  subroutine SetGeometryPointers ( G, iD )

    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iD

    select type ( C  =>  G % Geometry_C % Chart )
    class is ( Chart_BH_Form )

    call C % SetFieldPointer &
           ( G % Value ( :, G % EDGE_I_U ( iD ) ), Edge_I_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % WIDTH_U ( iD ) ),  Width_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % CENTER_U ( iD ) ), Center_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % AREA_I_D ( iD ) ), Area_I_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % VOLUME ), Volume_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % METRIC_F_DD_11 ), Metric_F_DD_11_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % METRIC_F_DD_22 ), Metric_F_DD_22_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % METRIC_F_DD_33 ), Metric_F_DD_33_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % METRIC_F_UU_11 ), Metric_F_UU_11_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % METRIC_F_UU_22 ), Metric_F_UU_22_3D )
    call C % SetFieldPointer &
           ( G % Value ( :, G % METRIC_F_UU_33 ), Metric_F_UU_33_3D )

    end select !-- C

  end subroutine SetGeometryPointers


end program Chart_BH__Form_Test
