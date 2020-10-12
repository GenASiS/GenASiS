module LaplacianMultipole_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: LaplacianMultipoleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0!, &
    character ( LDF ) :: &
      Type = '', &
      Name = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type LaplacianMultipoleTemplate

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

    ! call L % SetParameters ( A, MaxDegree, nEquations )
    ! call L % AllocateRectangularCoordinates ( )
    ! call L % AllocateSolidHarmonics ( )
    ! call L % AllocateMoments ( )

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( L )

    class ( LaplacianMultipoleTemplate ), intent ( inout ) :: &
      L

    ! if ( allocated ( L % ReductionMoments ) ) &
    !   deallocate ( L % ReductionMoments )

    ! if ( allocated ( L % RadialEdges ) ) &
    !   deallocate ( L % RadialEdges )
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


end module LaplacianMultipole_Template
