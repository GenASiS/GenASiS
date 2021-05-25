module InteractionsExamples_ASC__Form

  !-- InteractionsExamples_AtlasSingleChart_Form

  use GenASiS
  use Interactions_C__Form
  use Interactions_MWV_1__Form
  use Interactions_MWV_2__Form
  use Interactions_MWV_3__Form
  use InteractionsExamples_CSL__Form

  implicit none
  private

  type, public, extends ( Interactions_ASC_Template ) :: &
    InteractionsExamples_ASC_Form
  contains
    procedure, public, pass :: &
      Set_C_Grey
    procedure, private, pass :: &
      Set_MWV_1_Grey
    procedure, private, pass :: &
      Set_MWV_2_Grey
    procedure, private, pass :: &
      Set_MWV_3_Grey
    generic, public :: &
      Set_MWV_Grey => Set_MWV_1_Grey, Set_MWV_2_Grey, Set_MWV_3_Grey
    final :: &
      Finalize
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type InteractionsExamples_ASC_Form

contains


  subroutine Set_C_Grey ( IA, FA, OpacityAbsorption )

    class ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
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
      class is ( Interactions_C_Form )
        call I % Set ( Fluid = F, OpacityAbsorption = OpacityAbsorption )
      end select !-- I
    case ( 'SPECTRAL' )
      !-- Set at the Bundle level
    end select !-- MomentsType
    nullify ( F, I )

  end subroutine Set_C_Grey


  subroutine Set_MWV_1_Grey ( IA, FA, SpecificOpacity )

    class ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( InteractionsTemplate ), pointer :: &
      I

    select case ( trim ( IA % MomentsType ) )
    case ( 'GREY' )
      F => FA % Fluid_P_I ( )
      I => IA % Interactions ( )
      select type ( I )
      class is ( Interactions_MWV_1_Form )
        call I % Set ( Fluid = F, SpecificOpacity = SpecificOpacity )
      end select !-- I
    case ( 'SPECTRAL' )
      !-- Set at the Bundle level
    end select !-- MomentsType
    nullify ( F, I )

  end subroutine Set_MWV_1_Grey


  subroutine Set_MWV_2_Grey ( IA, FA, SpecificOpacity, EnergyMax )

    class ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      EnergyMax

    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( InteractionsTemplate ), pointer :: &
      I

    select case ( trim ( IA % MomentsType ) )
    case ( 'GREY' )
      F => FA % Fluid_P_I ( )
      I => IA % Interactions ( )
      select type ( I )
      class is ( Interactions_MWV_2_Form )
        call I % Set ( Fluid = F, SpecificOpacity = SpecificOpacity, &
                       EnergyMax = EnergyMax )
      end select !-- I
    case ( 'SPECTRAL' )
      !-- Set at the Bundle level
    end select !-- MomentsType
    nullify ( F, I )

  end subroutine Set_MWV_2_Grey


  subroutine Set_MWV_3_Grey &
               ( IA, FA, SpecificOpacity, EnergyMax, TemperatureScale )

    class ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      EnergyMax, &
      TemperatureScale

    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( InteractionsTemplate ), pointer :: &
      I

    select case ( trim ( IA % MomentsType ) )
    case ( 'GREY' )
      F => FA % Fluid_P_I ( )
      I => IA % Interactions ( )
      select type ( I )
      class is ( Interactions_MWV_3_Form )
        call I % Set ( Fluid = F, SpecificOpacity = SpecificOpacity, &
                       EnergyMax = EnergyMax, &
                       TemperatureScale = TemperatureScale )
      end select !-- I
    case ( 'SPECTRAL' )
      !-- Set at the Bundle level
    end select !-- MomentsType
    nullify ( F, I )

  end subroutine Set_MWV_3_Grey


  impure elemental subroutine Finalize ( IA )

    type ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
      IA

    call IA % FinalizeTemplate_I_ASC ( )

  end subroutine Finalize


  subroutine SetType ( IA )

    class ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
      IA

    if ( IA % Type == '' ) &
      IA % Type = 'an InteractionsExamples_ASC'

  end subroutine SetType


  subroutine AllocateField ( IA )

    class ( InteractionsExamples_ASC_Form ), intent ( inout ) :: &
      IA

    allocate ( InteractionsExamples_CSL_Form :: IA % Chart )

  end subroutine AllocateField


end module InteractionsExamples_ASC__Form
