module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use Mathematics
  use Fluid_ASC__Form
  use WoosleyHeger_07_Header_Form

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: WoosleyHeger_07_A_Form
    type ( WoosleyHeger_07_HeaderForm ), allocatable :: &
      Header
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type WoosleyHeger_07_A_Form

contains


  subroutine Initialize ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    allocate ( WH % Header )
    associate ( WHH => WH % Header )

    !-- PositionSpace

    call WHH % InitializePositionSpace ( WH )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: WH % Current_ASC )
    select type ( FA => WH % Current_ASC )
    class is ( Fluid_ASC_Form )

    call WHH % InitializeFluid ( FA )

    !-- Cleanup

    end select !-- FA
    end associate !-- WHH

  end subroutine Initialize


  impure elemental subroutine Finalize ( WH )

    type ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

    if ( allocated ( WH % Header ) ) &
      deallocate ( WH % Header )

    call WH % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


end module WoosleyHeger_07_A__Form
