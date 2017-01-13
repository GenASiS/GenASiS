module Interactions_CSL__Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_C__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Interactions_CSL_Form
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
  end type Interactions_CSL_Form

contains


  subroutine Initialize &
               ( IC, C, nValues, NameOutputOption, EquilibriumDensityOption, &
                 EffectiveOpacityOption, TransportOpacityOption )

    class ( Interactions_CSL_Form ), intent ( inout ) :: &
      IC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    real ( KDR ), intent ( in ), optional :: &
      EquilibriumDensityOption, &
      EffectiveOpacityOption, &
      TransportOpacityOption

    if ( IC % Type == '' ) &
      IC % Type = 'an Interactions_CSL'

    if ( IC % InteractionsType == '' ) &
      IC % InteractionsType = 'CONSTANT'    

    if ( present ( EquilibriumDensityOption ) ) &
      IC % EquilibriumDensity  =  EquilibriumDensityOption
    if ( present ( EffectiveOpacityOption ) ) &
      IC % EffectiveOpacity  =  EffectiveOpacityOption
    if ( present ( TransportOpacityOption ) ) &
      IC % TransportOpacity  =  TransportOpacityOption

    call IC % InitializeTemplate_CSL &
           ( C, nValues, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  subroutine SetField ( FC, NameOption )

    class ( Interactions_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate ( FC % FieldOutput )

    select case ( trim ( FC % InteractionsType ) )
    case ( 'CONSTANT' )
      allocate ( Interactions_C_Form :: FC % Field )
      select type ( I => FC % Field )
      type is ( Interactions_C_Form )
        call I % Initialize &
               ( FC % EquilibriumDensity, FC % EffectiveOpacity, &
                 FC % TransportOpacity, FC % nValues, NameOption = NameOption )
        call I % SetOutput ( FC % FieldOutput )
      end select !-- RM
    case default
      call Show ( 'InteractionsType not recognized', CONSOLE % ERROR )
      call Show ( FC % InteractionsType, 'InteractionsType', &
                  CONSOLE % ERROR )
      call Show ( 'Interactions_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- InteractionsType

  end subroutine SetField


  impure elemental subroutine Finalize ( IC )

    type ( Interactions_CSL_Form ), intent ( inout ) :: &
      IC

    call IC % FinalizeTemplate ( )

  end subroutine Finalize


end module Interactions_CSL__Form
