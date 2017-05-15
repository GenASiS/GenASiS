module Sources_F_CSL__Form

  !-- Sources_Fluid_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Sources_F__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Sources_F_CSL_Form
    type ( MeasuredValueForm ) :: &
      TimeUnit
    class ( Field_CSL_Template ), pointer :: &
      Fluid_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Sources_F_CSL_Form

contains


  subroutine Initialize &
               ( SFC, Fluid_CSL, NameShort, TimeUnit, nValues, &
                 IgnorabilityOption )

    class ( Sources_F_CSL_Form ), intent ( inout ) :: &
      SFC
    class ( Field_CSL_Template ), intent ( in ), target :: &
      Fluid_CSL
    character ( * ), intent ( in ) :: &
      NameShort
    type ( MeasuredValueForm ) :: &
      TimeUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( SFC % Type == '' ) &
      SFC % Type = 'a Sources_F_CSL'

    SFC % TimeUnit = TimeUnit

    SFC % Fluid_CSL => Fluid_CSL

    call SFC % InitializeTemplate_CSL &
           ( Fluid_CSL % Chart, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SFC )

    type ( Sources_F_CSL_Form ), intent ( inout ) :: &
      SFC

    nullify ( SFC % Fluid_CSL )

    call SFC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( Sources_F_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    allocate ( Sources_F_Form :: FC % Field )
    select type ( SF => FC % Field )
    class is ( Sources_F_Form )
    select type ( F => FC % Fluid_CSL % Field )
    class is ( Fluid_D_Form )
      call SF % Initialize ( F, FC % TimeUnit, NameOption = FC % NameShort )
      call SF % SetOutput ( FC % FieldOutput )
    end select !-- F
    end select !-- SF

  end subroutine SetField


end module Sources_F_CSL__Form
