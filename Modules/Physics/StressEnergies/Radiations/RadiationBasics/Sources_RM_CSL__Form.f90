module Sources_RM_CSL__Form

  !-- Sources_RadiationMoments_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Sources_RM_CSL_Form
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      EnergyUnit
    class ( Field_CSL_Template ), pointer :: &
      RadiationMoments_CSL => null ( )
  contains
    procedure, public, pass :: &
     Initialize
   final :: &
     Finalize
   procedure, private, pass :: &
     SetField
  end type Sources_RM_CSL_Form

contains


  subroutine Initialize &
               ( SRMC, RadiationMoments_CSL, NameShort, UsePinnedMemory, &
                 TimeUnit, EnergyUnit, nValues, IgnorabilityOption )

    class ( Sources_RM_CSL_Form ), intent ( inout ) :: &
      SRMC
    class ( Field_CSL_Template ), intent ( in ), target :: &
      RadiationMoments_CSL
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ) :: &
      UsePinnedMemory
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      EnergyUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( SRMC % Type == '' ) &
      SRMC % Type = 'a Sources_RM_CSL'

    SRMC % TimeUnit   = TimeUnit
    SRMC % EnergyUnit = EnergyUnit

    SRMC % RadiationMoments_CSL => RadiationMoments_CSL

    call SRMC % InitializeTemplate_CSL &
           ( RadiationMoments_CSL % Chart, NameShort, UsePinnedMemory, &
             nValues, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SRMC )

    type ( Sources_RM_CSL_Form ), intent ( inout ) :: &
      SRMC

    nullify ( SRMC % RadiationMoments_CSL )

    call SRMC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( Sources_RM_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    allocate ( Sources_RM_Form :: FC % Field )
    select type ( SRM => FC % Field )
    class is ( Sources_RM_Form )
    select type ( RM => FC % RadiationMoments_CSL % Field )
    class is ( RadiationMomentsForm )
      call SRM % Initialize &
             ( RM, FC % TimeUnit, FC % EnergyUnit, &
               PinnedOption = FC % UsePinnedMemory, &
               NameOption = FC % NameShort )
      call SRM % SetOutput ( FC % FieldOutput )
    end select !-- RM
    end select !-- SRM

  end subroutine SetField


end module Sources_RM_CSL__Form
