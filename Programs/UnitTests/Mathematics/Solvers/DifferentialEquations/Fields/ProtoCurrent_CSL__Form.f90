module ProtoCurrent_CSL__Form

  use Basics
  use Manifolds
  use Sources_C__Form
  use ProtoCurrent_Form
  use ProtoCurrentSources_CSL__Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: ProtoCurrent_CSL_Form
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit    
    class ( ProtoCurrentSources_CSL_Form ), pointer :: &
      Sources_CSL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      ProtoCurrent
    procedure, public, pass :: &
      SetSources
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type ProtoCurrent_CSL_Form

contains


  subroutine Initialize ( PCC, C, NameShort, nValues, IgnorabilityOption )

    class ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      PCC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( PCC % Type == '' ) &
      PCC % Type = 'a ProtoCurrent_CSL'

    call PCC % InitializeTemplate_CSL &
           ( C, NameShort, UsePinnedMemory = .true., &
             nValues = nValues, IgnorabilityOption = IgnorabilityOption )

  end subroutine Initialize


  function ProtoCurrent ( PCC ) result ( PC )

    class ( ProtoCurrent_CSL_Form ), intent ( in ), target :: &
      PCC
    class ( ProtoCurrentForm ), pointer :: &
      PC
      
    class ( StorageForm ), pointer :: &
      Field
    
    Field => PCC % Field
    select type ( Field )
    class is ( ProtoCurrentForm )
    PC => Field
    end select !-- Field

  end function ProtoCurrent


  subroutine SetSources ( PCC, PCSC )

    class ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      PCC
    class ( ProtoCurrentSources_CSL_Form ), intent ( in ), target :: &
      PCSC

    class ( ProtoCurrentForm ), pointer :: &
      PC

    PCC % Sources_CSL => PCSC

    PC => PCC % ProtoCurrent ( )
    select type ( PCS => PCSC % Field )
    class is ( Sources_C_Form )
      call PC % SetSources ( PCS )
    end select !-- FF

    nullify ( PC )

  end subroutine SetSources


  impure elemental subroutine Finalize ( PCC )

    type ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      PCC

    nullify ( PCC % Sources_CSL )

    call PCC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( ProtoCurrentForm :: FC % Field )
    allocate ( FC % FieldOutput )

    select type ( PC => FC % Field )
    class is ( ProtoCurrentForm )
      call PC % Initialize &
             ( FC % Velocity_U_Unit, FC % nValues, &
               NameOption = FC % NameShort )
      call PC % SetPrimitiveConserved ( )
      call PC % SetOutput ( FC % FieldOutput )
    end select !-- F

  end subroutine SetField


end module ProtoCurrent_CSL__Form
