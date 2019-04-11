module WoosleyHeger_07_D__Form

  !-- WoosleyHeger_07_Deleptonization__Form

  use GenASiS
  use WoosleyHeger_07_A__Form

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_A_Form ) :: WoosleyHeger_07_D_Form
  contains
    procedure, private, pass :: &
      Initialize_WH
  end type WoosleyHeger_07_D_Form

    private :: &
      ApplySources

contains


  subroutine Initialize_WH ( WH, Name )

    class ( WoosleyHeger_07_D_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    if ( WH % Type == '' ) &
      WH % Type = 'a WoosleyHeger_07_D'

    call WH % WoosleyHeger_07_A_Form % Initialize ( Name )

    select type ( I => WH % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( S => I % Step )
    class is ( Step_RK2_C_ASC_Form )
      S % ApplySources % Pointer => ApplySources
    end select !-- S

    end select !-- I

  end subroutine Initialize_WH


  subroutine ApplySources &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( StorageForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call ApplyCurvilinear_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    call ApplyGravity_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    call ApplyDeleptonization_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

  end subroutine ApplySources


end module WoosleyHeger_07_D__Form
