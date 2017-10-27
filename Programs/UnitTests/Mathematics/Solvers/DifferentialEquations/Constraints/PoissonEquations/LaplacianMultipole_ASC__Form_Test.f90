module LaplacianMultipole_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_ASC__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'LaplacianMultipole_ASC__Form_Test'
    
  type, public :: LaplacianMultipole_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      Stream
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
    type ( LaplacianMultipole_ASC_Form ), allocatable :: &
      Laplacian
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type LaplacianMultipole_ASC__Form_Test_Form

contains


  subroutine Initialize ( LMFT, Name )

    class ( LaplacianMultipole_ASC__Form_Test_Form ) :: &
      LMFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iA, iM, iEll, &  !-- iAngular, iOrder, iDegree
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal, &
      nCellsCore, &
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
    character ( LDL ), dimension ( : ), allocatable :: &
      R_C_Name, I_C_Name, &
      R_S_Name, I_S_Name
    type ( VariableGroupForm ) :: &
      Source, &
      SolidHarmonics_RC, SolidHarmonics_IC, &
      SolidHarmonics_RS, SolidHarmonics_IS
    class ( GeometryFlatForm ), pointer :: &
      G

    !-- Atlas

    allocate ( LMFT % Atlas )
    associate ( A => LMFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )

    CoordinateSystem = 'SPHERICAL'

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

    allocate ( LMFT % Geometry )
    associate ( GA => LMFT % Geometry )  !-- GeometryAtlas
    call GA % Initialize ( A )
    call A % SetGeometry ( GA )
    G => A % Geometry ( )
    end associate !-- GA

    !-- Laplacian

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( LMFT % Laplacian )
    associate ( L => LMFT % Laplacian )
    call L % Initialize ( A, MaxDegree )

    call Show ( L % RadialEdge, 'RadialEdge' )

    !-- Source

    call Source % Initialize &
           ( [ G % nValues, 1 ], &
               NameOption = 'Source', &
               VariableOption = [ 'HomogeneousDensity' ], &
               ClearOption = .true. )

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

    !-- Multipole names

    call L % NameSolidHarmonics ( R_C_Name, I_C_Name, R_S_Name, I_S_Name )
    call Show ( R_C_Name, '>>> R_C_Name' )
    call Show ( I_C_Name, '>>> I_C_Name' )
    if ( L % MaxOrder > 0 ) then
      call Show ( R_S_Name, '>>> R_S_Name' )
      call Show ( I_S_Name, '>>> I_S_Name' )
    end if

    !-- Compute solid harmonics

    call SolidHarmonics_RC % Initialize &
           ( [ G % nValues, size ( R_C_Name ) ], &
             NameOption = 'SolidHarmonics_RC', VariableOption = R_C_Name, &
             ClearOption = .true. )
    call SolidHarmonics_IC % Initialize &
           ( [ G % nValues, size ( I_C_Name ) ], &
             NameOption = 'SolidHarmonics_IC', VariableOption = I_C_Name, &
             ClearOption = .true. )
    if ( L % MaxOrder > 0 ) then
      call SolidHarmonics_RS % Initialize &
             ( [ G % nValues, size ( R_S_Name ) ], &
               NameOption = 'SolidHarmonics_RS', VariableOption = R_S_Name, &
               ClearOption = .true. )
      call SolidHarmonics_IS % Initialize &
             ( [ G % nValues, size ( I_S_Name ) ], &
               NameOption = 'SolidHarmonics_IS', VariableOption = I_S_Name, &
               ClearOption = .true. )
    end if

    !$OMP parallel do private ( iC )
    do iC = 1, G % nValues
      if ( .not. C % IsProperCell ( iC ) ) &
        cycle
      call L % ComputeSolidHarmonics &
             ( C % CoordinateSystem, &
               G % Value ( iC, G % CENTER ( 1 ) : G % CENTER ( 3 ) ), &
               C % nDimensions )
      SolidHarmonics_RC % Value ( iC, : )  =  L % SolidHarmonic_RC
      SolidHarmonics_IC % Value ( iC, : )  =  L % SolidHarmonic_IC
      if ( L % MaxOrder > 0 ) then
        SolidHarmonics_RS % Value ( iC, : )  =  L % SolidHarmonic_RS
        SolidHarmonics_IS % Value ( iC, : )  =  L % SolidHarmonic_IS
      end if
    end do
    !$OMP end parallel do

    !-- Compute moments

    call L % ComputeMoments ( Source )

    iA = 0
    do iM = 0, L % MaxOrder
      do iEll = iM, L % MaxDegree
        iA = iA + 1
        call Show ( iEll, 'iEll' )
        call Show ( iM, 'iM' )
        call Show ( L % MRC ( iA, : ), 'Moment_RC' )
        call Show ( L % MIC ( iA, : ), 'Moment_IC' )
        if ( iM > 0 ) then
          call Show ( L % MRS ( iA, : ), 'Moment_RS' )
          call Show ( L % MIS ( iA, : ), 'Moment_IS' )
        end if
      end do
    end do

    !-- Write

    allocate ( LMFT % Stream )
    call LMFT % Stream % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )    
    associate ( GIS => LMFT % Stream )

    call A % OpenStream ( GIS, 'Stream', iStream = 1 )
    call C % AddFieldImage ( Source, iStream = 1 )
    call C % AddFieldImage ( SolidHarmonics_RC, iStream = 1 )
    call C % AddFieldImage ( SolidHarmonics_IC, iStream = 1 )
    if ( L % MaxOrder > 0 ) then
      call C % AddFieldImage ( SolidHarmonics_RS, iStream = 1 )
      call C % AddFieldImage ( SolidHarmonics_IS, iStream = 1 )
    end if

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- GIS

    !-- Cleanup

    end associate !-- L
    end select !-- C
    end associate !-- A

    nullify ( G )

  end subroutine Initialize


  subroutine Finalize ( LMFT )

    type ( LaplacianMultipole_ASC__Form_Test_Form ) :: &
      LMFT

    if ( allocated ( LMFT % Laplacian ) ) &
      deallocate ( LMFT % Laplacian )
    if ( allocated ( LMFT % Geometry ) ) &
      deallocate ( LMFT % Geometry )
    if ( allocated ( LMFT % Atlas ) ) &
      deallocate ( LMFT % Atlas )
    if ( allocated ( LMFT % Stream ) ) &
      deallocate ( LMFT % Stream )

  end subroutine Finalize


end module LaplacianMultipole_ASC__Form_Test__Form


program LaplacianMultipole_ASC__Form_Test

  use Basics
  use LaplacianMultipole_ASC__Form_Test__Form

  implicit none

  type ( LaplacianMultipole_ASC__Form_Test_Form ), allocatable :: &
    LMFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( LMFT )
  call LMFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( LMFT )

  deallocate ( PROGRAM_HEADER )

end program LaplacianMultipole_ASC__Form_Test

