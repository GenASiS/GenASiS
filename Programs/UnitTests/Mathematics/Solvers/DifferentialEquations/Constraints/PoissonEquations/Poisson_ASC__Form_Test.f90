module Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Poisson_ASC__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Poisson_ASC__Form_Test'

  type, public :: Poisson_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      Stream
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
    type ( Storage_ASC_Form ), allocatable :: &
      Source, &
      Solution
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Poisson_ASC__Form_Test_Form

contains


  subroutine Initialize ( PFT, Name )

    class ( Poisson_ASC__Form_Test_Form ) :: &
      PFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      ! iC, &  !-- iCell
      ! iA, iM, iEll, &  !-- iAngular, iOrder, iDegree
      ! iR, & !-- iRadial
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal, &
      nCellsCore, &
      nEquations, &
      MaxDegree
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      RadiusMax, &
      RadiusCore, &
      RadiusDensity, &
      Density
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing
    type ( VariableGroupForm ), pointer :: &
      Source
    class ( GeometryFlatForm ), pointer :: &
      G

    !-- Atlas

    allocate ( PFT % Atlas )
    associate ( A => PFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )

    CoordinateSystem = 'SPHERICAL'

    call A % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
    if ( A % nDimensions > 1 ) &
      call A % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    if ( A % nDimensions > 2 ) &
      call A % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 3 )

    RadiusMax = 10.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'RadiusMax' )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [   0.0_KDR, 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    RadiusCore = RadiusMax / 8.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusCore, 'RadiusCore' )

    nCellsCore = 32  !-- Number of central cells with equal spacing
    call PROGRAM_HEADER % GetParameter ( nCellsCore, 'nCellsCore' )

    call Show ( 'Mesh core parameters' )
    call Show ( RadiusCore, 'RadiusCore' )
    call Show ( nCellsCore, 'nCellsCore' )
    call Show ( RadiusCore / nCellsCore, 'CellWidthCore' )

    nCellsRadial = 3 * nCellsCore  !-- Aiming for roughly R_max = 10
    call PROGRAM_HEADER % GetParameter ( nCellsRadial, 'nCellsRadial' )

        nCellsPolar = 3 * nCellsCore
    nCellsAzimuthal = 2 * nCellsPolar
 
    nCells = [ nCellsRadial, 1, 1 ]
    if ( A % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( A % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  RadiusCore

    call A % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nEqualOption = nCellsCore )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )


    !-- Geometry

    allocate ( PFT % Geometry )
    associate ( GA => PFT % Geometry )  !-- GeometryAtlas
    call GA % Initialize ( A )
    call A % SetGeometry ( GA )
    G => A % Geometry ( )
    end associate !-- GA


    !-- Poisson

    nEquations = 1

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( PFT % Poisson )
    associate ( P => PFT % Poisson )
    call P % Initialize &
           ( A, SolverType = 'MULTIPOLE', MaxDegreeOption = MaxDegree, &
             nEquationsOption = nEquations )


    !-- Source

    allocate ( PFT % Source )
    associate ( SA => PFT % Source )

    call SA % Initialize &
           ( A, 'Source', nEquations, &
             VariableOption = [ 'HomogeneousDensity' ], &
             WriteOption = .true. )
    Source => SA % Storage ( )

    RadiusDensity = RadiusMax / 10.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    associate &
      ( R => G % Value ( :, G % CENTER ( 1 ) ), &
        S => Source % Value ( :, 1 ) )
    where ( R < RadiusDensity )
      S = Density
    end where
    end associate !-- R, etc.

    end associate !-- SA


    !-- Write

    allocate ( PFT % Stream )
    call PFT % Stream % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )    
    associate ( GIS => PFT % Stream )

    call A % OpenStream ( GIS, 'Stream', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- GIS


    !-- Cleanup

    end associate !-- P
    end select !-- C
    end associate !-- A

    nullify ( G, Source )

  end subroutine Initialize


  subroutine Finalize ( PFT )

    type ( Poisson_ASC__Form_Test_Form ) :: &
      PFT

    if ( allocated ( PFT % Stream ) ) &
      deallocate ( PFT % Stream )
    if ( allocated ( PFT % Solution ) ) &
      deallocate ( PFT % Solution )
    if ( allocated ( PFT % Source ) ) &
      deallocate ( PFT % Source )
    if ( allocated ( PFT % Poisson ) ) &
      deallocate ( PFT % Poisson )
    if ( allocated ( PFT % Geometry ) ) &
      deallocate ( PFT % Geometry )
    if ( allocated ( PFT % Atlas ) ) &
      deallocate ( PFT % Atlas )

  end subroutine Finalize


end module Poisson_ASC__Form_Test__Form



program Poisson_ASC__Form_Test

  use Basics
  use Poisson_ASC__Form_Test__Form

  implicit none

  type ( Poisson_ASC__Form_Test_Form ), allocatable :: &
    PFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( PFT )
  call PFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( PFT )

  deallocate ( PROGRAM_HEADER )

end program Poisson_ASC__Form_Test
