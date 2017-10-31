module Poisson_Template

  use Basics
  use Manifolds

  implicit none
  private

  type, public, abstract :: PoissonTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Type = '', &
      Name = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type PoissonTemplate


contains


  subroutine InitializeTemplate ( P, A )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    P % IGNORABILITY = A % IGNORABILITY

    if ( P % Type == '' ) &
      P % Type = 'a Poisson' 

    P % Name = 'Poisson_' // trim ( A % Name )

    call Show ( 'Initializing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( P )

    class ( PoissonTemplate ), intent ( inout ) :: &
      P

    if ( P % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( P % Type ), P % IGNORABILITY )
    call Show ( P % Name, 'Name', P % IGNORABILITY )
    
  end subroutine FinalizeTemplate


end module Poisson_Template
