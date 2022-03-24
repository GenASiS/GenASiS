module Series_F_CC__Form

  !-- Series_Fluid_CentralCore__Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public, extends ( Series_CS_Form ) :: Series_F_CC_Form
  contains
    final :: &
      Finalize
  end type Series_F_CC_Form


contains


  impure elemental subroutine Finalize ( S )

    type ( Series_F_CC_Form ), intent ( inout ) :: &
      S

    ! if ( allocated ( S % Change ) ) &
    !   deallocate ( S % Change )
    ! if ( allocated ( S % Total ) ) &
    !   deallocate ( S % Total )
    ! if ( allocated ( S % Boundary ) ) &
    !   deallocate ( S % Boundary )
    ! if ( allocated ( S % Interior ) ) &
    !   deallocate ( S % Interior )

    ! nullify ( S % TallyChange )
    ! nullify ( S % TallyTotal )
    ! nullify ( S % TallyBoundary )
    ! nullify ( S % TallyInterior )

  end subroutine Finalize


end module Series_F_CC__Form
