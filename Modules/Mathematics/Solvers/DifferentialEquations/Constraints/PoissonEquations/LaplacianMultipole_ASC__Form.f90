module LaplacianMultipole_ASC__Form

  !-- LaplacianMultipole_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, extends ( LaplacianMultipoleTemplate ) :: &
    LaplacianMultipole_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type LaplacianMultipole_ASC_Form

contains


  subroutine Initialize ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole_ASC' 

    call L % InitializeTemplate ( A, MaxDegree, nEquations )

  end subroutine Initialize


  impure elemental subroutine Finalize ( L )

    type ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    ! nullify ( L % Chart )

    ! if ( allocated ( L % SolidHarmonics ) ) &
    !   deallocate ( L % SolidHarmonics )
    ! if ( allocated ( L % RectangularCoordinates ) ) &
    !   deallocate ( L % RectangularCoordinates )

    ! nullify ( L % Source )
    ! nullify ( L % SolidHarmonic_IS )
    ! nullify ( L % SolidHarmonic_RS )
    ! nullify ( L % SolidHarmonic_IC )
    ! nullify ( L % SolidHarmonic_RC )

    ! nullify ( L % Volume )
    ! nullify ( L % RadiusSquared )
    ! nullify ( L % Radius )
    ! nullify ( L % Rectangular_Z )
    ! nullify ( L % Rectangular_Y )
    ! nullify ( L % Rectangular_X )

    call L % FinalizeTemplate ( )

  end subroutine Finalize


end module LaplacianMultipole_ASC__Form
