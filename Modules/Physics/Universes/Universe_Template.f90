module Universe_Template

  use Basics
  use Mathematics
 
  implicit none
  private

  type, public, abstract :: UniverseTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    class ( IntegratorTemplate ), allocatable :: &
      Integrator
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type UniverseTemplate

contains


  subroutine InitializeTemplate ( U, Name )

    class ( UniverseTemplate ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      Name

    U % IGNORABILITY = CONSOLE % INFO_1

    if ( U % Type == '' ) &
      U % Type = 'a Universe' 

    U % Name = Name

    call Show ( 'Initializing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( U )

    class ( UniverseTemplate ), intent ( inout ) :: &
      U

    if ( allocated ( U % Integrator ) ) &
      deallocate ( U % Integrator )

    nullify ( U % Communicator )

    if ( U % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

  end subroutine FinalizeTemplate


end module Universe_Template