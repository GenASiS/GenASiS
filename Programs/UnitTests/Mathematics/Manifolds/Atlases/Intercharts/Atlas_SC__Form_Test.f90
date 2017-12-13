module Atlas_SC__Form_Test__Form

  use Basics
  use Charts
  use GeometryFlat_ASC__Form
  use Atlas_SC__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Atlas_SC__Form_Test'

  type, public :: Atlas_SC_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GIS_Base
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      GeometryFiber
    type ( Atlas_SC_Form ), allocatable :: &
      AtlasBase, &
      AtlasFiber
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Atlas_SC_Form_Test_Form

contains

  
  subroutine Initialize ( AFT, Name )

    class ( Atlas_SC_Form_Test_Form ), intent ( inout ), target :: &
      AFT
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
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
      ( AFT % AtlasBase, &
        AFT % GIS_Base )
    associate &
      ( A   => AFT % AtlasBase, &
        GIS => AFT % GIS_Base )

    call A % Initialize &
           ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
             iDimensionalityOption = 1 )
    call A % CreateChart ( )  
    call A % SetGeometry ( )

    call GIS % Initialize &
           ( PROGRAM_HEADER % Name, CommunicatorOption = A % Communicator )
    call A % OpenStream ( GIS, '1', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- A, etc.

    !-- Fiber Atlas

    allocate &
      ( AFT % AtlasFiber, &
        AFT % GeometryFiber )
    associate &
      ( A => AFT % AtlasFiber, &
        G => AFT % GeometryFiber )

    call A % Initialize ( 'Fiber', iDimensionalityOption = 2 )
    call A % SetBoundaryConditionsFace ( [ 'REFLECTING', 'REFLECTING' ], 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = 5.0_KDR * UNIT % MEGA_ELECTRON_VOLT

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEGA_ELECTRON_VOLT

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    call A % CreateChart &
           ( SpacingOption = Spacing, CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, nGhostLayersOption = [ 0, 0, 0 ] )

    call G % InitializeFlat ( A )
    call A % SetGeometry ( GeometryOption = G )

    associate ( GIS_Base => AFT % GIS_Base )

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
    call A % OpenStream ( GIS_Fiber, '1', iStream = 1 )

    call GIS_Fiber % Open &
           ( GIS_Fiber % ACCESS_CREATE, SeriesOption = .false. )
    call A % Write ( iStream = 1 )
    call GIS_Fiber % Close ( )

    call A % CloseStreams ( )
    deallocate ( GIS_Fiber )
    
    end associate !-- GIS_Base
    end associate !-- A, etc.

  end subroutine Initialize


  subroutine Finalize ( AFT )

    type ( Atlas_SC_Form_Test_Form ), intent ( inout ) :: &
      AFT

    deallocate ( AFT % AtlasFiber )
    deallocate ( AFT % AtlasBase )
    deallocate ( AFT % GeometryFiber )
    deallocate ( AFT % GIS_Base )

  end subroutine Finalize


end module Atlas_SC__Form_Test__Form



program Atlas_SC__Form_Test

  use Basics
  use Atlas_SC__Form_Test__Form
  
  implicit none

  type ( Atlas_SC_Form_Test_Form ), allocatable :: &
    AFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '3D_1D' )

  allocate ( AFT )
  call AFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( AFT )

  deallocate ( PROGRAM_HEADER )

end program Atlas_SC__Form_Test
