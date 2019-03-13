module Universe_Template

  use Basics
  use Mathematics
  use StressEnergies
 
  implicit none
  private

  type, public, extends ( UniverseHeaderForm ), abstract :: UniverseTemplate
    type ( StressEnergyUnitsForm ) :: &
      Units
    class ( IntegratorTemplate ), allocatable :: &
      Integrator
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      Evolve
    procedure, public, pass :: &
      FinalizeTemplate
  end type UniverseTemplate

contains


  subroutine InitializeTemplate ( U, Name )

    class ( UniverseTemplate ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      Name

    call U % InitializeHeader ( Name )

  end subroutine InitializeTemplate


  subroutine Evolve ( U )

    class ( UniverseTemplate ), intent ( inout ) :: &
      U

    call U % Integrator % Evolve ( )

  end subroutine Evolve


  impure elemental subroutine FinalizeTemplate ( U )

    class ( UniverseTemplate ), intent ( inout ) :: &
      U

    if ( allocated ( U % Integrator ) ) &
      deallocate ( U % Integrator )

  end subroutine FinalizeTemplate


end module Universe_Template
