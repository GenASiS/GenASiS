module FluidFeatures_ASC__Form

  !-- FluidFeatures_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use FluidFeatures_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: FluidFeatures_ASC_Form
    real ( KDR ) :: &
      ShockThreshold, &
      TrivialDensity
    character ( LDF ) :: &
      FluidType = ''
    class ( Field_ASC_Template ), pointer :: &
      Fluid_ASC => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type FluidFeatures_ASC_Form

contains


  subroutine Initialize &
               ( FFA, Fluid_ASC, FluidType, NameShortOption, &
                 ShockThresholdOption, TrivialDensityOption, &
                 IgnorabilityOption )

    class ( FluidFeatures_ASC_Form ), intent ( inout ) :: &
      FFA
    class ( Field_ASC_Template ), intent ( in ), target :: &
      Fluid_ASC
    character ( * ), intent ( in ) :: &
      FluidType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    real ( KDR ), intent ( in ), optional :: &
      ShockThresholdOption, &
      TrivialDensityOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort

    if ( FFA % Type == '' ) &
      FFA % Type = 'a FluidFeatures_ASC'
    FFA % FluidType = FluidType

    FFA % ShockThreshold = 0.2_KDR
    if ( present ( ShockThresholdOption ) ) &
      FFA % ShockThreshold = ShockThresholdOption
    call PROGRAM_HEADER % GetParameter &
           ( FFA % ShockThreshold, 'ShockThreshold' ) 

    FFA % TrivialDensity = 0.0_KDR
    if ( present ( TrivialDensityOption ) ) &
      FFA % TrivialDensity = TrivialDensityOption
    call PROGRAM_HEADER % GetParameter &
           ( FFA % TrivialDensity, 'TrivialDensity' ) 

    FFA % Fluid_ASC => Fluid_ASC

    NameShort = 'FluidFeatures'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call FFA % InitializeTemplate_ASC &
           ( Fluid_ASC % Atlas, NameShort, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FFA )

    type ( FluidFeatures_ASC_Form ), intent ( inout ) :: &
      FFA

    nullify ( FFA % Fluid_ASC )

    call FFA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( FluidFeatures_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    select type ( FC => FA % Fluid_ASC % Chart )
    class is ( Field_CSL_Template )

    allocate ( FluidFeatures_CSL_Form :: FA % Chart )

    select type ( FFC => FA % Chart )
    class is ( FluidFeatures_CSL_Form )
      call FFC % Initialize &
             ( FC, FA % NameShort, FA % FluidType, FA % ShockThreshold, &
               FA % TrivialDensity, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FFC

    call A % AddField ( FA )

    end select !-- FC
    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module FluidFeatures_ASC__Form
