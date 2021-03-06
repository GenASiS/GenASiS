module Chart_SLX__Form_Test__Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLL__Form
  use Chart_SLD__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Chart_SLX__Form_Test'

  type, public :: Chart_SLX__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GIS_Base
    type ( AtlasHeaderForm ), allocatable :: &
      AtlasHeaderBase, &
      AtlasHeaderFiber
    type ( GeometryFlat_CSL_Form ), allocatable :: &
      GeometryBase, &
      GeometryFiber
    type ( Chart_SLL_Form ), allocatable :: &
      ChartFiber
    type ( Chart_SLD_Form ), allocatable :: &
      ChartBase
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Chart_SLX__Form_Test_Form

contains

  
  subroutine Initialize ( CFT )

    class ( Chart_SLX__Form_Test_Form ), intent ( inout ), target :: &
      CFT

    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      G_3D
    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic
    character ( LDN + 1 ) :: &
      FileNumber, &
      BlockNumber
    character ( LDF ) :: &
      BundleDirectory, &
      BundleName
    character ( LDL ), dimension ( 3 ) :: &
      Spacing
    type ( GridImageStreamForm ), allocatable :: &
      GIS_Fiber

    !-- Base Atlas 

    allocate &
      ( CFT % AtlasHeaderBase, &
        CFT % GeometryBase, &
        CFT % ChartBase, &
        CFT % GIS_Base )
    associate &
      ( A => CFT % AtlasHeaderBase, &
        GC => CFT % GeometryBase, &
        C => CFT % ChartBase, &
        GIS => CFT % GIS_Base )

    call A % Initialize &
           ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
             iDimensionalityOption = 1 )

    IsPeriodic = .true.

    call C % Initialize ( A, IsPeriodic, iChart = 1 )

    associate ( nValues => C % nProperCells + C % nGhostCells )
    call GC % InitializeFlat &
           ( C, 'Geometry', UsePinnedMemory = .false., nValues = nValues )
    call C % AddField ( GC )
    C % iFieldGeometry = C % nFields
    call C % SetGeometry ( GC )
    end associate !-- nValues

    call A % Show ( )
    call C % Show ( )

    call GIS % Initialize &
           ( PROGRAM_HEADER % Name, CommunicatorOption = A % Communicator )
    call C % OpenStream ( GIS, '1', iStream = 1 )

    select type ( G  =>  GC % Field )
    class is ( GeometryFlatForm )

    call C % SetVariablePointer ( G % Value ( :, G % CENTER_U ( 1 ) ), G_3D )
    call Show ( G_3D ( :, 3, 3 ), 'Center_1 reference' ) !-- Indexed for 3D

    call GIS % Open ( GIS % ACCESS_CREATE )
    call C % Write ( iStream = 1 )
    call GIS % Close ( )

    associate ( nCB => C % nCellsBrick )
    G_3D ( 1 : nCB ( 1 ), 1 : nCB ( 2 ), 1 : nCB ( 3 ) )  =  0.0_KDR
    end associate !-- nCB
    call Show ( G_3D ( :, 3, 3 ), 'Center_1 proper cleared before read' )

    call GIS % Open ( GIS % ACCESS_READ, NumberOption = GIS % Number )
    call C % Read ( iStream = 1 )
    call GIS % Close ( )

    call Show ( G_3D ( :, 3, 3 ), 'Center_1 after read' )

    end select !-- G
    end associate !-- A, etc.

    !-- Fiber Atlas

    allocate &
      ( CFT % AtlasHeaderFiber, &
        CFT % GeometryFiber, &
        CFT % ChartFiber )
    associate &
      ( A => CFT % AtlasHeaderFiber, &
        GC => CFT % GeometryFiber, &
        C => CFT % ChartFiber )

    call A % Initialize ( 'Fiber', iDimensionalityOption = 2 )
    call A % SetBoundaryConditionsFace ( [ 'REFLECTING', 'REFLECTING' ], 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = 5.0_KDR * UNIT % MEGA_ELECTRON_VOLT

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEGA_ELECTRON_VOLT

    IsPeriodic = .false.

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    call C % Initialize &
           ( A, IsPeriodic, iChart = 1, SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, nGhostLayersOption = [ 0, 0, 0 ] )

    associate ( nValues => C % nProperCells + C % nGhostCells )
    call GC % InitializeFlat &
           ( C, 'Geometry', UsePinnedMemory = .false., nValues = nValues )
    call C % AddField ( GC ) 
    C % iFieldGeometry = C % nFields
    call C % SetGeometry ( GC )
    end associate !-- nValues

    call A % Show ( )
    call C % Show ( )

    associate ( GIS_Base => CFT % GIS_Base )

    write ( FileNumber, fmt = '(a1,i7.7)' ) '_', &
      GIS_Base % Number    
    write ( BlockNumber, fmt = '(a1,i7.7)' ) '_', &
      GIS_Base % Communicator % Rank
    BundleDirectory &
      = trim ( GIS_Base % WorkingDirectory ) // trim ( GIS_Base % Name ) &
        // '_Bundle_MeshBlock'// trim ( BlockNumber ) // '/' 
    BundleName &
      = trim ( GIS_Base % Name ) // '_Bundle_MeshBlock' &
        // trim ( BlockNumber ) // trim ( FileNumber )
    call Show ( BundleDirectory, 'BundleDirectory' )
    call Show ( BundleName, 'BundleName' )

    allocate ( GIS_Fiber )
    call GIS_Fiber % Initialize &
           ( BundleName, WorkingDirectoryOption = BundleDirectory )
    call C % OpenStream ( GIS_Fiber, '1', iStream = 1 )

    call GIS_Fiber % Open &
           ( GIS_Fiber % ACCESS_CREATE, SeriesOption = .false. )
    call C % Write ( iStream = 1 )
    call GIS_Fiber % Close ( )

    call C % CloseStreams ( )
    deallocate ( GIS_Fiber )
    
    end associate !-- GIS_Base
    end associate !-- A, etc.

  end subroutine Initialize


  subroutine Finalize ( CFT )

    type ( Chart_SLX__Form_Test_Form ), intent ( inout ) :: &
      CFT

    deallocate ( CFT % ChartFiber )
    deallocate ( CFT % AtlasHeaderFiber )

    deallocate ( CFT % GIS_Base )
    deallocate ( CFT % ChartBase )
    deallocate ( CFT % GeometryBase )
    deallocate ( CFT % AtlasHeaderBase )

  end subroutine Finalize


end module Chart_SLX__Form_Test__Form



program Chart_SLX__Form_Test

  use Basics
  use Chart_SLX__Form_Test__Form
  
  implicit none

  type ( Chart_SLX__Form_Test_Form ), allocatable :: &
    CFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '3D_1D' )
    
  allocate ( CFT )
  call CFT % Initialize ( )
  deallocate ( CFT )

  deallocate ( PROGRAM_HEADER )

end program Chart_SLX__Form_Test
