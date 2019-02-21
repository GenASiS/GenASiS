module Interactions_G_ASC__Form

  !-- Interactions_Generic_AtlasSingleChart_Form

  use GenASiS
  use Interactions_G__Form
  use Interactions_G_CSL__Form

  implicit none
  private

  type, public, extends ( Interactions_ASC_Template ) :: &
    Interactions_G_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Set => Set_G_ASC_G
    procedure, public, pass :: &
      AllocateField
  end type Interactions_G_ASC_Form

contains


  subroutine Initialize &
               ( IA, A, InteractionsType, MomentsType, NameShortOption, &
                 LengthUnitOption, EnergyDensityUnitOption, &
                 TemperatureUnitOption, IgnorabilityOption )

    class ( Interactions_G_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      InteractionsType, &
      MomentsType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      EnergyDensityUnitOption, &
      TemperatureUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( IA % Type == '' ) &
      IA % Type = 'an Interactions_G_ASC'

    call IA % InitializeTemplate_I_ASC &
           ( A, InteractionsType, MomentsType, NameShortOption, &
             LengthUnitOption, EnergyDensityUnitOption, TemperatureUnitOption, &
             IgnorabilityOption )

  end subroutine Initialize


  subroutine Set_G_ASC_G ( IA, FA, OpacityAbsorption )

    class ( Interactions_G_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      OpacityAbsorption

    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( InteractionsTemplate ), pointer :: &
      I

    select case ( trim ( IA % MomentsType ) )
    case ( 'GREY' )
      F => FA % Fluid_P_I ( )
      I => IA % Interactions ( )
      select type ( I )
      type is ( Interactions_G_Form )
        call I % Set ( Fluid = F, OpacityAbsorption = OpacityAbsorption )
      end select !-- I
    case ( 'SPECTRAL' )
      !-- Set at the Bundle level
    end select !-- MomentsType
    nullify ( F, I )

  end subroutine Set_G_ASC_G


  subroutine AllocateField ( IA )

    class ( Interactions_G_ASC_Form ), intent ( inout ) :: &
      IA

    allocate ( Interactions_G_CSL_Form :: IA % Chart )

  end subroutine AllocateField


end module Interactions_G_ASC__Form
