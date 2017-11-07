module LaplacianMultipole_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_ASC__Form
  use CreateProportionalChart_Command
  use SetHomogeneousSphere_Command

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'LaplacianMultipole_ASC__Form_Test'
    
  type, public :: LaplacianMultipole_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      Stream
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
      nEquations, &
      MaxDegree
    real ( KDR ) :: &
      RadiusDensity, &
      Density, &
      Radius, &
      M_RC, M_IC, &
      M_RS, M_IS
    character ( 3 ) :: &
      Label_Ell, &
      Label_M
    character ( LDL ), dimension ( 1 ) :: &
      Variable
    character ( LDL ), dimension ( : ), allocatable :: &
      Suffix_Ell_M
    type ( VariableGroupForm ), pointer :: &
      Source, &
      SolidHarmonics_RC, SolidHarmonics_IC, &
      SolidHarmonics_RS, SolidHarmonics_IS, &
      Solution
    class ( GeometryFlatForm ), pointer :: &
      G


    !-- Atlas

    allocate ( LMFT % Atlas )
    associate ( A => LMFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )

    call CreateProportionalChart ( A )

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )


    !-- Laplacian

    nEquations = 1

    MaxDegree = 2
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( LMFT % Laplacian )
    associate ( L => LMFT % Laplacian )
    call L % Initialize ( A, MaxDegree, nEquations )

    call Show ( L % RadialEdge, 'RadialEdge', CONSOLE % INFO_2 )


    !-- Homogeneous sphere

    Variable = [ 'HomogeneousSphere' ]

    allocate ( LMFT % Source )
    associate ( SA => LMFT % Source )
    call SA % Initialize &
           ( A, 'Source', nEquations, &
             VariableOption = Variable, &
             WriteOption = .true. )
    Source => SA % Storage ( )

    allocate ( LMFT % Reference )
    associate ( RA => LMFT % Reference )
    call RA % Initialize &
           ( A, 'Reference', nEquations, &
             VariableOption = Variable, &
             WriteOption = .true. )

    RadiusDensity = C % MaxCoordinate ( 1 ) / 10.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    call SetHomogeneousSphere &
           ( SA, RA, A, Density, RadiusDensity, iVariable = 1 )

    end associate !-- RA
    end associate !-- SA


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
               C % nDimensions, Radius, iR )
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
        call Show ( L % M_RC ( iA, :, 1 ), 'Moment_RC', CONSOLE % INFO_3 )
        call Show ( L % M_IC ( iA, :, 1 ), 'Moment_IC', CONSOLE % INFO_3 )
        if ( iM > 0 ) then
          call Show ( L % M_RS ( iA, :, 1 ), 'Moment_RS', CONSOLE % INFO_3 )
          call Show ( L % M_IS ( iA, :, 1 ), 'Moment_IS', CONSOLE % INFO_3 )
        end if
      end do
    end do


    !-- Solution

    allocate ( LMFT % Solution )
    associate ( SA => LMFT % Solution )
    call SA % Initialize &
           ( A, 'Solution', L % nAngularMomentCells + 1, &
             VariableOption &
               = trim ( Variable ( 1 ) ) // [ Suffix_Ell_M, '' ], &
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
               C % nDimensions, Radius, iR )

      associate ( Phi => Solution % Value ( iC, L % nAngularMomentCells + 1 ) )
      do iA = 1, L % nAngularMomentCells
        associate ( Phi_Ell_M => Solution % Value ( iC, iA ) )

        M_RC  =  0.5_KDR  *  L % M_RC ( iA, iR, 1 )
        if ( iR > 1 ) &
          M_RC  =  M_RC  +  0.5_KDR  *  L % M_RC ( iA, iR - 1, 1 )

        M_IC  =  0.5_KDR  *  L % M_IC ( iA, iR, 1 )
        if ( iR < L % nRadialCells ) &
          M_IC  =  M_IC  +  0.5_KDR  *  L % M_IC ( iA, iR + 1, 1 )

        associate &
          ( R_C  =>  L % SolidHarmonic_RC ( iA ), &
            I_C  =>  L % SolidHarmonic_IC ( iA ) )
        Phi_Ell_M  =  L % Delta ( iA ) * ( M_RC * I_C  +  M_IC * R_C )  
        end associate !-- R_C, etc.
        
        if ( L % MaxOrder > 0 ) then

          M_RS  =  0.5_KDR  *  L % M_RS ( iA, iR, 1 )
          if ( iR > 1 ) &
            M_RS  =  M_RS  +  0.5_KDR  *  L % M_RS ( iA, iR - 1, 1 )

          M_IS  =  0.5_KDR  *  L % M_IS ( iA, iR, 1 )
          if ( iR < L % nRadialCells ) &
            M_IS  =  M_IS  +  0.5_KDR  *  L % M_IS ( iA, iR + 1, 1 )

          associate &
            ( R_S  =>  L % SolidHarmonic_RS ( iA ), &
              I_S  =>  L % SolidHarmonic_IS ( iA ) )
          Phi_Ell_M  =  Phi_Ell_M &
                        +  L % Delta ( iA ) * ( M_RS * I_S  +  M_IS * R_S )  
          end associate !-- R_C, etc.
        
        end if

        Phi_Ell_M  =  - Phi_Ell_M  /  ( 4 * CONSTANT % PI )
        Phi        =    Phi  +  Phi_Ell_M

        end associate !-- Phi_L_M
      end do !-- iA
      end associate !-- Phi

    end do !-- iC
    !$OMP end parallel do


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
              SolidHarmonics_RS, SolidHarmonics_IS, Solution )

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

