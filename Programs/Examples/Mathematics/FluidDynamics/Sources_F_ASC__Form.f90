module Sources_F_ASC__Form

  !-- Sources_Fluid_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use Sources_F_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Sources_F_ASC_Form
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
  end type Sources_F_ASC_Form

contains


  subroutine Initialize &
               ( SFA, Fluid_ASC, NameShortOption, TimeUnitOption, &
                 IgnorabilityOption )

    class ( Sources_F_ASC_Form ), intent ( inout ) :: &
      SFA
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

    if ( SFA % Type == '' ) &
      SFA % Type = 'a Sources_F_ASC'

    if ( present ( TimeUnitOption ) ) &
      SFA % TimeUnit = TimeUnitOption

    SFA % Fluid_ASC => Fluid_ASC

    NameShort = 'Sources_F'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call SFA % InitializeTemplate_ASC &
           ( Fluid_ASC % Atlas, NameShort, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SFA )

    type ( Sources_F_ASC_Form ), intent ( inout ) :: &
      SFA

    nullify ( SFA % Fluid_ASC )

    call SFA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Sources_F_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    select type ( FC => FA % Fluid_ASC % Chart )
    class is ( Field_CSL_Template )

    allocate ( Sources_F_CSL_Form :: FA % Chart )

    select type ( SFC => FA % Chart )
    class is ( Sources_F_CSL_Form )
      call SFC % Initialize &
             ( FC, FA % NameShort, FA % TimeUnit, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- SFC

    call A % AddField ( FA )

    end select !-- FC
    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Sources_F_ASC__Form
