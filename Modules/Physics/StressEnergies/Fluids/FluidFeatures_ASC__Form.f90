module FluidFeatures_ASC__Form

  !-- FluidFeatures_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use FluidFeatures_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: FluidFeatures_ASC_Form
    real ( KDR ) :: &
      ShockThreshold
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
               ( FFA, Fluid_ASC, FluidType, RiemannSolverType, &
                 NameShortOption, UsePinnedMemoryOption, &
                 ShockThresholdOption, IgnorabilityOption )

    class ( FluidFeatures_ASC_Form ), intent ( inout ) :: &
      FFA
    class ( Field_ASC_Template ), intent ( in ), target :: &
      Fluid_ASC
    character ( * ), intent ( in ) :: &
      FluidType, &
      RiemannSolverType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    real ( KDR ), intent ( in ), optional :: &
      ShockThresholdOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort

    if ( FFA % Type == '' ) &
      FFA % Type = 'a FluidFeatures_ASC'
    FFA % FluidType = FluidType

    FFA % ShockThreshold = 0.1_KDR
    if ( present ( ShockThresholdOption ) ) &
      FFA % ShockThreshold = ShockThresholdOption
    call PROGRAM_HEADER % GetParameter &
           ( FFA % ShockThreshold, 'ShockThreshold' ) 

    FFA % Fluid_ASC => Fluid_ASC

    NameShort = 'FluidFeatures'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call FFA % InitializeTemplate_ASC &
           ( Fluid_ASC % Atlas, NameShort, UsePinnedMemoryOption, &
             IgnorabilityOption )

    call Show ( FFA % ShockThreshold, 'ShockThreshold', FFA % IGNORABILITY )

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
             ( FC, FA % NameShort, FA % FluidType, FA % UsePinnedMemory, &
               FA % ShockThreshold, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FFC

    call A % AddField ( FA )

    end select !-- FC
    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module FluidFeatures_ASC__Form
