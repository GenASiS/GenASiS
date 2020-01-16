module Sources_RM_ASC__Form

  !-- Sources_RadiationMoments_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use Sources_RM_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Sources_RM_ASC_Form
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      EnergyUnit
    logical ( KDL ) :: &
      SuppressWrite
    class ( Field_ASC_Template ), pointer :: &
      RadiationMoments_ASC => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Sources_RM_ASC_Form

contains


  subroutine Initialize &
               ( SRMA, RadiationMoments_ASC, NameShortOption, TimeUnitOption, &
                 EnergyUnitOption, IgnorabilityOption, SuppressWriteOption )

    class ( Sources_RM_ASC_Form ), intent ( inout ) :: &
      SRMA
    class ( Field_ASC_Template ), intent ( in ), target :: &
      RadiationMoments_ASC
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption, &
      EnergyUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( in ), optional :: &
      SuppressWriteOption

    character ( LDL ) :: &
      NameShort

    if ( SRMA % Type == '' ) &
      SRMA % Type = 'a Sources_RM_ASC'

    if ( present ( TimeUnitOption ) ) &
      SRMA % TimeUnit = TimeUnitOption
    if ( present ( EnergyUnitOption ) ) &
      SRMA % EnergyUnit = EnergyUnitOption

    SRMA % RadiationMoments_ASC => RadiationMoments_ASC

    NameShort = 'Sources_RM'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    SRMA % SuppressWrite = .false.
    if ( present ( SuppressWriteOption ) ) &
      SRMA % SuppressWrite = SuppressWriteOption

    call SRMA % InitializeTemplate_ASC &
           ( RadiationMoments_ASC % Atlas, NameShort, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SRMA )

    type ( Sources_RM_ASC_Form ), intent ( inout ) :: &
      SRMA

    nullify ( SRMA % RadiationMoments_ASC )

    call SRMA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Sources_RM_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    select type ( RMC => FA % RadiationMoments_ASC % Chart )
    class is ( Field_CSL_Template )

    allocate ( Sources_RM_CSL_Form :: FA % Chart )

    select type ( SRMC => FA % Chart )
    class is ( Sources_RM_CSL_Form )
      call SRMC % Initialize &
             ( RMC, FA % NameShort, FA % TimeUnit, FA % EnergyUnit, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- SRMC

    if ( .not. FA % SuppressWrite ) &
      call A % AddField ( FA )

    end select !-- RMC
    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Sources_RM_ASC__Form
