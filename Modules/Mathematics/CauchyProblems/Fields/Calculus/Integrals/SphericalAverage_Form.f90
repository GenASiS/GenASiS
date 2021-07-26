module SphericalAverage_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries

  implicit none
  private

  type, public, extends ( FieldSetForm ) :: SphericalAverageForm
    class ( FieldSetForm ), pointer :: &
      Source
  contains
!    procedure, public, nopass :: &
!      Compute
    final :: &
      Finalize
  end type SphericalAverageForm


contains


  subroutine Initialize ( SA, Source )

    class ( SphericalAverageForm ), intent ( inout ) :: &
      SA
    class ( FieldSetForm ), intent ( in ), target :: &
      Source
    
    SA % Source  =>  Source

  end subroutine Initialize


!   subroutine Compute ( SA, I, G, IgnorabilityOption )

!     class ( FieldSetForm ), intent ( inout ) :: &
!       SA  !-- SphericalAverage
!     class ( FieldSetForm ), intent ( in ) :: &
!       I   !-- Integrand
!   end subroutine Compute


  impure elemental subroutine Finalize ( SA )

    type ( SphericalAverageForm ), intent ( inout ) :: &
      SA

  end subroutine Finalize


end module SphericalAverage_Form
