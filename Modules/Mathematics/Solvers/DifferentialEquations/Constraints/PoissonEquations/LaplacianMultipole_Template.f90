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
    logical ( KDL ) :: &
      UseDevice = .false., &
      ReductionUseDevice = .false.
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( StorageForm ), allocatable :: &
!        Moments, &
!      MyMoments, &
      RadialEdges
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure, private, pass :: &
      SetParameters
    procedure ( SPA ), private, pass, deferred :: &
      SetParametersAtlas
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

    L % IGNORABILITY  =  A % IGNORABILITY

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole' 

    L % Name = 'Laplacian_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )

    call L % SetParameters ( A, MaxDegree, nEquations )
    ! call L % AllocateRectangularCoordinates ( )
    ! call L % AllocateSolidHarmonics ( )
    ! call L % AllocateMoments ( )

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    ! if ( allocated ( L % ReductionMoments ) ) &
    !   deallocate ( L % ReductionMoments )

    if ( allocated ( L % RadialEdges ) ) &
      deallocate ( L % RadialEdges )
    ! if ( allocated ( L % MyMoments ) ) &
    !   deallocate ( L % MyMoments )
    ! if ( allocated ( L % Moments ) ) &
    !   deallocate ( L % Moments )

    ! nullify ( L % MyMoment_IS )
    ! nullify ( L % MyMoment_RS )
    ! nullify ( L % MyMoment_IC )
    ! nullify ( L % MyMoment_RC )
    ! nullify ( L % Moment_IS )
    ! nullify ( L % Moment_RS )
    ! nullify ( L % Moment_IC )
    ! nullify ( L % Moment_RC )

    if ( L % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( L % Type ), L % IGNORABILITY )
    call Show ( L % Name, 'Name', L % IGNORABILITY )
    
  end subroutine FinalizeTemplate


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
    call Show ( L % ReductionUseDevice, 'ReductionUseDevice', L % IGNORABILITY )

  end subroutine SetParameters


end module LaplacianMultipole_Template
