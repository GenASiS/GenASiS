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
    type ( Storage_ASC_Form ), allocatable :: &
      Source, &
      SolidHarmonics_RC, SolidHarmonics_IC, &
      SolidHarmonics_RS, SolidHarmonics_IS, &
      Solution, &
      Reference
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
      iR, & !-- iRadial
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
      Density, &
      MRC, MIC, &
      MRS, MIS
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    character ( 3 ) :: &
      Label_Ell, &
      Label_M
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing
    character ( LDL ), dimension ( : ), allocatable :: &
      Suffix_Ell_M
    type ( VariableGroupForm ), pointer :: &
      Source, &
      SolidHarmonics_RC, SolidHarmonics_IC, &
      SolidHarmonics_RS, SolidHarmonics_IS, &
      Solution, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G


    !-- Atlas

    allocate ( LMFT % Atlas )
    associate ( A => LMFT % Atlas )
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

    allocate ( LMFT % Geometry )
    associate ( GA => LMFT % Geometry )  !-- GeometryAtlas
    call GA % Initialize ( A )
    call A % SetGeometry ( GA )
    G => A % Geometry ( )
    end associate !-- GA


    !-- Laplacian

    nEquations = 1

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( LMFT % Laplacian )
    associate ( L => LMFT % Laplacian )
    call L % Initialize ( A, MaxDegree, nEquations )

    call Show ( L % RadialEdge, 'RadialEdge', CONSOLE % INFO_2 )


    !-- Source

    allocate ( LMFT % Source )
    associate ( SA => LMFT % Source )
    call SA % Initialize &
           ( A, 'Source', nEquations, &
             VariableOption = [ 'HomogeneousDensity' ], &
             WriteOption = .true. )
    Source => SA % Storage ( )
    end associate !-- SA

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


    !-- Multipole name suffixes

    allocate ( Suffix_Ell_M ( L % nAngularMomentCells ) )
    iA = 0
    do iM = 0, L % MaxOrder
      do iEll = iM, L % MaxDegree
        iA = iA + 1
        write ( Label_Ell, fmt = '(i3.3)' ) iEll
        write ( Label_M,   fmt = '(i3.3)' ) iM
        Suffix_Ell_M ( iA )  =  '_' // Label_Ell // '_' // Label_M
      end do
    end do


    !-- Compute solid harmonics

    call Show ( 'RegularCos' // Suffix_Ell_M, 'RC' )
    call Show ( 'IrregularCos' // Suffix_Ell_M, 'IC' )

    allocate ( LMFT % SolidHarmonics_RC, LMFT % SolidHarmonics_IC )
    associate &
      ( SH_RC => LMFT % SolidHarmonics_RC, &
        SH_IC => LMFT % SolidHarmonics_IC )
    call SH_RC % Initialize &
           ( A, 'SolidHarmonics_RC', L % nAngularMomentCells, &
             VariableOption = 'RegularCos' // Suffix_Ell_M, &
             WriteOption = .true. )
    call SH_IC % Initialize &
           ( A, 'SolidHarmonics_IC', L % nAngularMomentCells, &
             VariableOption = 'IrregularCos' // Suffix_Ell_M, &
             WriteOption = .true. )
    SolidHarmonics_RC => SH_RC % Storage ( )
    SolidHarmonics_IC => SH_IC % Storage ( )
    end associate !-- SH_RC, etc.
    
    if ( L % MaxOrder > 0 ) then

      call Show ( 'RegularSin' // Suffix_Ell_M, 'RS' )
      call Show ( 'IrregularSin' // Suffix_Ell_M, 'IS' )

      allocate ( LMFT % SolidHarmonics_RS, LMFT % SolidHarmonics_IS )
      associate &
        ( SH_RS => LMFT % SolidHarmonics_RS, &
          SH_IS => LMFT % SolidHarmonics_IS )
      call SH_RS % Initialize &
             ( A, 'SolidHarmonics_RS', L % nAngularMomentCells, &
               VariableOption = 'RegularSin' // Suffix_Ell_M, &
               WriteOption = .true. )
      call SH_IS % Initialize &
             ( A, 'SolidHarmonics_IS', L % nAngularMomentCells, &
               VariableOption = 'IrregularSin' // Suffix_Ell_M, &
               WriteOption = .true. )
      SolidHarmonics_RS => SH_RS % Storage ( )
      SolidHarmonics_IS => SH_IS % Storage ( )
      end associate !-- SH_RS, etc.
    
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
        call Show ( iEll, 'iEll', CONSOLE % INFO_3 )
        call Show ( iM, 'iM', CONSOLE % INFO_3 )
        call Show ( L % MRC ( iA, :, 1 ), 'Moment_RC', CONSOLE % INFO_3 )
        call Show ( L % MIC ( iA, :, 1 ), 'Moment_IC', CONSOLE % INFO_3 )
        if ( iM > 0 ) then
          call Show ( L % MRS ( iA, :, 1 ), 'Moment_RS', CONSOLE % INFO_3 )
          call Show ( L % MIS ( iA, :, 1 ), 'Moment_IS', CONSOLE % INFO_3 )
        end if
      end do
    end do


    !-- Solution

    allocate ( LMFT % Solution )
    associate ( SA => LMFT % Solution )
    call SA % Initialize &
           ( A, 'Solution', L % nAngularMomentCells + 1, &
             VariableOption = 'Phi' // [ Suffix_Ell_M, '' ], &
             WriteOption = .true. )
    Solution => SA % Storage ( )
    end associate !-- SA

    !$OMP parallel do private ( iC )
    do iC = 1, G % nValues

      if ( .not. C % IsProperCell ( iC ) ) &
        cycle

      call L % ComputeSolidHarmonics &
             ( C % CoordinateSystem, &
               G % Value ( iC, G % CENTER ( 1 ) : G % CENTER ( 3 ) ), &
               C % nDimensions )

      call Search ( L % RadialEdge, G % Value ( iC, G % CENTER ( 1 ) ), iR )

      associate ( Phi => Solution % Value ( iC, L % nAngularMomentCells + 1 ) )
      do iA = 1, L % nAngularMomentCells
        associate ( Phi_Ell_M => Solution % Value ( iC, iA ) )

        MRC  =  0.5_KDR  *  L % MRC ( iA, iR, 1 )
        if ( iR > 1 ) &
          MRC  =  MRC  +  0.5_KDR  *  L % MRC ( iA, iR - 1, 1 )

        MIC  =  0.5_KDR  *  L % MIC ( iA, iR, 1 )
        if ( iR < L % nRadialCells ) &
          MIC  =  MIC  +  0.5_KDR  *  L % MIC ( iA, iR + 1, 1 )

        associate &
          ( R_C  =>  L % SolidHarmonic_RC ( iA ), &
            I_C  =>  L % SolidHarmonic_IC ( iA ) )
        Phi_Ell_M  =  L % Delta ( iA ) * ( MRC * I_C  +  MIC * R_C )  
        end associate !-- R_C, etc.
        
        if ( L % MaxOrder > 0 ) then

          MRS  =  0.5_KDR  *  L % MRS ( iA, iR, 1 )
          if ( iR > 1 ) &
            MRS  =  MRS  +  0.5_KDR  *  L % MRS ( iA, iR - 1, 1 )

          MIS  =  0.5_KDR  *  L % MIS ( iA, iR, 1 )
          if ( iR < L % nRadialCells ) &
            MIS  =  MIS  +  0.5_KDR  *  L % MIS ( iA, iR + 1, 1 )

          associate &
            ( R_S  =>  L % SolidHarmonic_RS ( iA ), &
              I_S  =>  L % SolidHarmonic_IS ( iA ) )
          Phi_Ell_M  =  Phi_Ell_M &
                        +  L % Delta ( iA ) * ( MRS * I_S  +  MIS * R_S )  
          end associate !-- R_C, etc.
        
        end if

        Phi_Ell_M  =  - Phi_Ell_M  /  ( 4 * CONSTANT % PI )
        Phi        =    Phi  +  Phi_Ell_M

        end associate !-- Phi_L_M
      end do !-- iA
      end associate !-- Phi

    end do !-- iC
    !$OMP end parallel do


    !-- Reference

    allocate ( LMFT % Reference )
    associate ( RA => LMFT % Reference )
    call RA % Initialize &
           ( A, 'Reference', nEquations, &
             VariableOption = [ 'Phi' ], &
             WriteOption = .true. )
    Reference => RA % Storage ( )
    end associate !-- RA

    associate &
      ( R    =>  G % Value ( :, G % CENTER ( 1 ) ), &
        Phi  =>  Reference % Value ( :, 1 ), &
        Pi   =>  CONSTANT % PI )
    where ( R < RadiusDensity )
      Phi  =  1.0_KDR / 6.0_KDR  *  Density  *  R ** 2  &
              -  1.0_KDR / 2.0_KDR  *  Density  *  RadiusDensity ** 2
    elsewhere
      Phi  =  - 1.0_KDR / 3.0_KDR  *  Density  *  RadiusDensity ** 3  /  R
    end where
    end associate !-- R


    !-- Write

    allocate ( LMFT % Stream )
    call LMFT % Stream % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )    
    associate ( GIS => LMFT % Stream )

    call A % OpenStream ( GIS, 'Stream', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- GIS

    !-- Cleanup

    end associate !-- L
    end select !-- C
    end associate !-- A

    nullify ( G, Source, SolidHarmonics_RC, SolidHarmonics_IC, &
              SolidHarmonics_RS, SolidHarmonics_IS, Solution, Reference )

  end subroutine Initialize


  subroutine Finalize ( LMFT )

    type ( LaplacianMultipole_ASC__Form_Test_Form ) :: &
      LMFT

    if ( allocated ( LMFT % Stream ) ) &
      deallocate ( LMFT % Stream )
    if ( allocated ( LMFT % Reference ) ) &
      deallocate ( LMFT % Reference )
    if ( allocated ( LMFT % Solution ) ) &
      deallocate ( LMFT % Solution )
    if ( allocated ( LMFT % SolidHarmonics_IS ) ) &
      deallocate ( LMFT % SolidHarmonics_IS )
    if ( allocated ( LMFT % SolidHarmonics_RS ) ) &
      deallocate ( LMFT % SolidHarmonics_RS )
    if ( allocated ( LMFT % SolidHarmonics_IC ) ) &
      deallocate ( LMFT % SolidHarmonics_IC )
    if ( allocated ( LMFT % SolidHarmonics_RC ) ) &
      deallocate ( LMFT % SolidHarmonics_RC )
    if ( allocated ( LMFT % Source ) ) &
      deallocate ( LMFT % Source )
    if ( allocated ( LMFT % Laplacian ) ) &
      deallocate ( LMFT % Laplacian )
    if ( allocated ( LMFT % Geometry ) ) &
      deallocate ( LMFT % Geometry )
    if ( allocated ( LMFT % Atlas ) ) &
      deallocate ( LMFT % Atlas )

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

