module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use Mathematics
  use Fluid_P_MHN__Form
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
    procedure, public, pass :: &
      SetWriteTimeInterval
    procedure, public, pass :: &
      PrepareCycle
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

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: WH % Step )
    select type ( S => WH % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplyGravity

    !-- Initial Conditions

    call WHH % SetFluid ( )

    !-- Integrator

    call WHH % SetIntegratorParameters ( )

    call WH % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = WHH % TimeUnit, &
             FinishTimeOption = WHH % FinishTime, &
             CourantFactorOption = WHH % CourantFactor, &
             nWriteOption = WHH % nWrite )

    !-- Cleanup

    end select !-- S
    end select !-- FA
    end associate !-- WHH

  end subroutine Initialize


  subroutine SetWriteTimeInterval ( I )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      I

    call I % Header % SetWriteTimeInterval ( )

  end subroutine SetWriteTimeInterval


  subroutine PrepareCycle ( I )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      I

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_MHN_Form ), pointer :: &
      F

    associate ( WHH => I % Header )

    associate ( FA => WHH % Fluid_ASC )
    F => FA % Fluid_P_MHN ( )

    select type ( PS => WHH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call WHH % ComputeGravity ( F, PS, G )

    end select !-- PS
    end associate !-- FA
    end associate !-- WHH
    nullify ( F, G )

  end subroutine PrepareCycle


  impure elemental subroutine Finalize ( WH )

    type ( WoosleyHeger_07_A_Form ), intent ( inout ) :: &
      WH

    if ( allocated ( WH % Header ) ) &
      deallocate ( WH % Header )

    call WH % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


end module WoosleyHeger_07_A__Form
