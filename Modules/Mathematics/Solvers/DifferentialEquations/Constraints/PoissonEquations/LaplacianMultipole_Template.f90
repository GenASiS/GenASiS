module LaplacianMultipole_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: LaplacianMultipoleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nRadialCells = 0, &
      nAngularMoments = 0, &
      nEquations = 0, &
      MaxDegree = 0, &  !-- Max L
      MaxOrder  = 0     !-- Max M
    integer ( KDI ) :: &
      REGULAR_COSINE    = 1, &
      IRREGULAR_COSINE  = 2, &
      REGULAR_SINE      = 3, &
      IRREGULAR_SINE    = 4
    integer ( KDI ) :: &
      iSolidHarmonics, &
      iSolidHarmonics_P1, &
      iSolidHarmonics_P2
    real ( KDR ), dimension ( 3 ) :: &
      Origin = 0.0_KDR
    logical ( KDL ) :: &
      UseDevice = .false.
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), pointer :: &
      SolidHarmonics    => null ( ), &
      SolidHarmonics_P1 => null ( ), &  !-- Previous_1
      SolidHarmonics_P2 => null ( ), &  !-- Previous_2
      SolidHarmonics_PD => null ( )     !-- PreviousDiagonal
    type ( StorageForm ), dimension ( : ), allocatable :: &
      SolidHarmonics_1D
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, private, pass :: &
      SetParameters
    procedure ( SPA ), private, pass, deferred :: &
      SetParametersAtlas
    procedure ( ASH ), private, pass, deferred :: &
      AllocateSolidHarmonics
    procedure, public, pass :: &
      FinalizeTemplate
  end type LaplacianMultipoleTemplate


  abstract interface

    subroutine SPA ( L, A )
      use Manifolds
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
      class ( AtlasHeaderForm ), intent ( in ), target :: &
        A
    end subroutine SPA

    subroutine ASH ( L )
      import LaplacianMultipoleTemplate
      implicit none
      class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
        L
    end subroutine ASH

  end interface


contains


  subroutine InitializeTemplate ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    L % IGNORABILITY = A % IGNORABILITY

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole' 

    L % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call L % SetParameters ( A, MaxDegree, nEquations )
    call L % AllocateSolidHarmonics ( )

call Show ( '>>> Aborting during development', CONSOLE % ERROR )
call Show ( 'LaplacianMultipole_Template', 'module', CONSOLE % ERROR )
call Show ( 'InitializeTemplate', 'subroutine', CONSOLE % ERROR )
call PROGRAM_HEADER % Abort ( )

  end subroutine InitializeTemplate


  subroutine SetParameters ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    integer ( KDI ) :: &
      iL, &
      iM

    L % MaxDegree  =  MaxDegree
    L % MaxOrder   =  MaxDegree

    call L % SetParametersAtlas ( A )

    associate &
      (  L  =>  L % MaxDegree, &
         M  =>  L % MaxOrder, &
        nA  =>  L % nAngularMoments )
    nA = 0
    do iM  =  0, M
      do iL  =  iM, L
        nA  =  nA + 1
      end do
    end do
    end associate !-- L, etc.

    L % nEquations  =  nEquations

    call Show ( L % MaxDegree, 'MaxDegree (l)', L % IGNORABILITY )
    call Show ( L % MaxOrder, 'MaxOrder (m)', L % IGNORABILITY )
    call Show ( L % nRadialCells, 'nRadialCells', L % IGNORABILITY )
    call Show ( L % nAngularMoments, 'nAngularMoments', L % IGNORABILITY )
    call Show ( L % nEquations, 'nEquations', L % IGNORABILITY )
    call Show ( L % UseDevice, 'UseDevice', L % IGNORABILITY )

  end subroutine SetParameters


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    if ( allocated ( L % SolidHarmonics_1D ) ) &
      deallocate ( L % SolidHarmonics_1D )

    nullify ( L % SolidHarmonics_PD )
    nullify ( L % SolidHarmonics_P2 )
    nullify ( L % SolidHarmonics_P1 )
    nullify ( L % SolidHarmonics )

    if ( L % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )
    
  end subroutine FinalizeTemplate


end module LaplacianMultipole_Template
