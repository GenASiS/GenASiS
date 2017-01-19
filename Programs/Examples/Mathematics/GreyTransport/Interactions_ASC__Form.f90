module Interactions_ASC__Form

  use Basics
  use Mathematics
  use Interactions_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Interactions_ASC_Form
    real ( KDR ) :: &  !-- parameters for CONSTANT Interactions
      EquilibriumDensity = 0.0_KDR, &
      EffectiveOpacity   = 0.0_KDR, &
      TransportOpacity   = 0.0_KDR
    character ( LDL ) :: &
      InteractionsType = ''
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Interactions_ASC_Form

contains


  subroutine Initialize &
               ( IA, A, InteractionsType, NameOutputOption, &
                 EquilibriumDensityOption, EffectiveOpacityOption, &
                 TransportOpacityOption )

    class ( Interactions_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      InteractionsType
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    real ( KDR ), intent ( in ), optional :: &
      EquilibriumDensityOption, &
      EffectiveOpacityOption, &
      TransportOpacityOption

    if ( IA % Type == '' ) &
      IA % Type = 'an Interactions_ASC'
    IA % InteractionsType = InteractionsType    

    if ( present ( EquilibriumDensityOption ) ) &
      IA % EquilibriumDensity  =  EquilibriumDensityOption
    if ( present ( EffectiveOpacityOption ) ) &
      IA % EffectiveOpacity  =  EffectiveOpacityOption
    if ( present ( TransportOpacityOption ) ) &
      IA % TransportOpacity  =  TransportOpacityOption

    call IA % InitializeTemplate_ASC &
           ( A, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( IA )

    type ( Interactions_ASC_Form ), intent ( inout ) :: &
      IA

    call IA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA, NameOutputOption )

    class ( Interactions_ASC_Form ), intent ( inout ) :: &
      FA
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Interactions_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Interactions_CSL_Form )
      select case ( trim ( FA % InteractionsType ) )
      case ( 'CONSTANT' )
        call FC % Initialize &
               ( C, FA % InteractionsType, nValues, &
                 NameOutputOption = NameOutputOption, &
                 EquilibriumDensityOption = FA % EquilibriumDensity, &
                 EffectiveOpacityOption   = FA % EffectiveOpacity, &
                 TransportOpacityOption   = FA % TransportOpacity )
      end select !-- InteractionsType
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Interactions_ASC__Form
