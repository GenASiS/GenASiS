module LaplacianMultipole_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: LaplacianMultipoleForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      MaxMoment
    character ( LDF ) :: &
      Name = ''
    class ( ChartTemplate ), pointer :: &
      Chart => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type LaplacianMultipoleForm

contains


  subroutine Initialize ( LM, Chart, MaxMoment )

    class ( LaplacianMultipoleForm ), intent ( inout ) :: &
      LM
    class ( ChartTemplate ), intent ( in ), target :: &
      Chart
    integer ( KDI ), intent ( in ) :: &
      MaxMoment

    LM % IGNORABILITY = Chart % IGNORABILITY
    LM % Name = 'Laplacian_' // trim ( Chart % Name )

    call Show ( 'Initializing a LaplacianMultipole', LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )
    call Show ( MaxMoment, 'MaxMoment', LM % IGNORABILITY )

    LM % MaxMoment = MaxMoment
    LM % Chart => Chart

  end subroutine Initialize


  impure elemental subroutine Finalize ( LM )

    type ( LaplacianMultipoleForm ), intent ( inout ) :: &
      LM

    nullify ( LM % Chart )

    if ( LM % Name == '' ) return

    call Show ( 'Finalizing a LaplacianMultipole', LM % IGNORABILITY )
    call Show ( LM % Name, 'Name', LM % IGNORABILITY )
    
  end subroutine Finalize


end module LaplacianMultipole_Form
