module FluidSources_CSL__Form

  !-- FluidSources_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use Fluid_D__Form
  use FluidSources_Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: FluidSources_CSL_Form
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
  end type FluidSources_CSL_Form

contains


  subroutine Initialize &
               ( FSC, Fluid_CSL, NameShort, TimeUnit, nValues, &
                 IgnorabilityOption )

    class ( FluidSources_CSL_Form ), intent ( inout ) :: &
      FSC
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

    if ( FSC % Type == '' ) &
      FSC % Type = 'a Fluid_CSL'

    FSC % TimeUnit = TimeUnit

    FSC % Fluid_CSL => Fluid_CSL

    call FSC % InitializeTemplate_CSL &
           ( Fluid_CSL % Chart, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FSC )

    type ( FluidSources_CSL_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % Fluid_CSL )

    call FSC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( FluidSources_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    allocate ( FluidSourcesForm :: FC % Field )
    select type ( FS => FC % Field )
    class is ( FluidSourcesForm )
    select type ( F => FC % Fluid_CSL % Field )
    class is ( Fluid_D_Form )
      call FS % Initialize &
             ( F, FC % TimeUnit, F % iaConserved, &
               NameOption = FC % NameShort )
      call FS % SetOutput ( FC % FieldOutput )
    end select !-- F
    end select !-- FS

  end subroutine SetField


end module FluidSources_CSL__Form
