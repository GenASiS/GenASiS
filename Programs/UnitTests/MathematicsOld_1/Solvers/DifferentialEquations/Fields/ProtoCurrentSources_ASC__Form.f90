module ProtoCurrentSources_ASC__Form

  !-- ProtoCurrentSources_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use ProtoCurrentSources_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: ProtoCurrentSources_ASC_Form
    type ( MeasuredValueForm ) :: &
      TimeUnit
    class ( Field_ASC_Template ), pointer :: &
      ProtoCurrent_ASC => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type ProtoCurrentSources_ASC_Form

contains


  subroutine Initialize &
               ( PCSA, ProtoCurrent_ASC, NameShortOption, TimeUnitOption, &
                 IgnorabilityOption )

    class ( ProtoCurrentSources_ASC_Form ), intent ( inout ) :: &
      PCSA
    class ( Field_ASC_Template ), intent ( in ), target :: &
      ProtoCurrent_ASC
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      NameShort

    if ( PCSA % Type == '' ) &
      PCSA % Type = 'a ProtoCurrentSources_ASC'

    if ( present ( TimeUnitOption ) ) &
      PCSA % TimeUnit = TimeUnitOption

    PCSA % ProtoCurrent_ASC => ProtoCurrent_ASC

    NameShort = 'ProtoCurrentSources'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call PCSA % InitializeTemplate_ASC &
           ( ProtoCurrent_ASC % Atlas, NameShort, &
             UsePinnedMemoryOption = .true., &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( PCSA )

    type ( ProtoCurrentSources_ASC_Form ), intent ( inout ) :: &
      PCSA

    nullify ( PCSA % ProtoCurrent_ASC )

    call PCSA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( ProtoCurrentSources_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    select type ( PCC => FA % ProtoCurrent_ASC % Chart )
    class is ( Field_CSL_Template )

    allocate ( ProtoCurrentSources_CSL_Form :: FA % Chart )

    select type ( PCSC => FA % Chart )
    class is ( ProtoCurrentSources_CSL_Form )
      call PCSC % Initialize &
             ( PCC, FA % NameShort, FA % TimeUnit, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- PCSC

    call A % AddField ( FA )

    end select !-- PCC
    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module ProtoCurrentSources_ASC__Form
