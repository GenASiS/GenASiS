module Poisson_ASC__Form

  !-- Poisson_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use Poisson_Template

  implicit none
  private

  type, public, extends ( PoissonTemplate ) :: Poisson_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Poisson_ASC_Form


contains


  subroutine Initialize ( P, A )

    class ( Poisson_ASC_Form ), intent ( inout ) :: &
      P
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A

    if ( P % Type == '' ) &
      P % Type = 'a Poisson_ASC' 

    call P % InitializeTemplate ( A )

  end subroutine Initialize


  impure elemental subroutine Finalize ( P )

    type ( Poisson_ASC_Form ), intent ( inout ) :: &
      P

    call P % FinalizeTemplate ( )

  end subroutine Finalize


end module Poisson_ASC__Form
