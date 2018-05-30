module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use GenASiS
  use WoosleyHeger_07__Template

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Template ) :: WoosleyHeger_07_A_Form
  contains
    procedure, public, pass :: &
      Initialize
  end type WoosleyHeger_07_A_Form

contains


  subroutine Initialize ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    if ( WH % Type == '' ) &
      WH % Type = 'a WoosleyHeger_07_A'

    call WH % InitializeTemplate ( Name )

    !-- Integrator

    allocate ( FluidCentralCoreForm :: WH % Integrator )
    select type ( FCC => WH % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'HEAVY_NUCLEUS', &
             GeometryType = 'NEWTONIAN', &
             ShockThresholdOption = 1.0_KDR, &
             nWriteOption = 30 )

    !-- Initial Conditions

    call WH % SetFluid ( )

    !-- Cleanup

    end select !-- FCC

  end subroutine Initialize


end module WoosleyHeger_07_A__Form
