module ProtoCurrentSources_CSL__Form

  !-- ProtoCurrentSources_ChartSingleLevel__Form

  use Basics
  use Manifolds
  use CurrentSources_Form
  use ProtoCurrent_Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: ProtoCurrentSources_CSL_Form
    type ( MeasuredValueForm ) :: &
      TimeUnit
    class ( Field_CSL_Template ), pointer :: &
      ProtoCurrent_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type ProtoCurrentSources_CSL_Form

contains


  subroutine Initialize &
               ( PCSC, ProtoCurrent_CSL, NameShort, TimeUnit, nValues, &
                 IgnorabilityOption )

    class ( ProtoCurrentSources_CSL_Form ), intent ( inout ) :: &
      PCSC
    class ( Field_CSL_Template ), intent ( in ), target :: &
      ProtoCurrent_CSL
    character ( * ), intent ( in ) :: &
      NameShort
    type ( MeasuredValueForm ) :: &
      TimeUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( PCSC % Type == '' ) &
      PCSC % Type = 'a ProtoCurrentSources_CSL'

    PCSC % TimeUnit = TimeUnit

    PCSC % ProtoCurrent_CSL => ProtoCurrent_CSL

    call PCSC % InitializeTemplate_CSL &
           ( ProtoCurrent_CSL % Chart, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( PCSC )

    type ( ProtoCurrentSources_CSL_Form ), intent ( inout ) :: &
      PCSC

    nullify ( PCSC % ProtoCurrent_CSL )

    call PCSC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( ProtoCurrentSources_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )

    allocate ( CurrentSourcesForm :: FC % Field )
    select type ( PCS => FC % Field )
    class is ( CurrentSourcesForm )
    select type ( PC => FC % ProtoCurrent_CSL % Field )
    class is ( ProtoCurrentForm )
      call PCS % Initialize &
             ( PC, FC % TimeUnit, PC % iaConserved, &
               NameOption = FC % NameShort )
      call PCS % SetOutput ( FC % FieldOutput )
    end select !-- PC
    end select !-- PCS

  end subroutine SetField


end module ProtoCurrentSources_CSL__Form
