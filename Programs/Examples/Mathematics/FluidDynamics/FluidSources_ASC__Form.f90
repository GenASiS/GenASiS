module FluidSources_ASC__Form

  !-- FluidSources_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use FluidSources_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: FluidSources_ASC_Form
    type ( MeasuredValueForm ) :: &
      TimeUnit
    class ( Field_ASC_Template ), pointer :: &
      Fluid_ASC => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type FluidSources_ASC_Form

contains


  subroutine Initialize &
               ( FSA, Fluid_ASC, NameShortOption, TimeUnitOption, &
                 IgnorabilityOption )

    class ( FluidSources_ASC_Form ), intent ( inout ) :: &
      FSA
    class ( Field_ASC_Template ), intent ( in ), target :: &
      Fluid_ASC
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      NameShort

    if ( FSA % Type == '' ) &
      FSA % Type = 'a FluidSources_ASC'

    if ( present ( TimeUnitOption ) ) &
      FSA % TimeUnit = TimeUnitOption

    FSA % Fluid_ASC => Fluid_ASC

    NameShort = 'FluidSources'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call FSA % InitializeTemplate_ASC &
           ( Fluid_ASC % Atlas, NameShort, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FSA )

    type ( FluidSources_ASC_Form ), intent ( inout ) :: &
      FSA

    nullify ( FSA % Fluid_ASC )

    call FSA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( FluidSources_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    select type ( FC => FA % Fluid_ASC % Chart )
    class is ( Field_CSL_Template )

    allocate ( FluidSources_CSL_Form :: FA % Chart )

    select type ( FSC => FA % Chart )
    class is ( FluidSources_CSL_Form )
      call FSC % Initialize &
             ( FC, FA % NameShort, FA % TimeUnit, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FSC

    call A % AddField ( FA )

    end select !-- FC
    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module FluidSources_ASC__Form
