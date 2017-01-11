module ProtoCurrent_CSL__Form

  use Basics
  use Manifolds
  use ProtoCurrent_Form

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: ProtoCurrent_CSL_Form
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit    
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      ProtoCurrent
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type ProtoCurrent_CSL_Form

contains


  subroutine Initialize ( PCC, C, nValues, NameOutputOption )

    class ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      PCC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    if ( PCC % Type == '' ) &
      PCC % Type = 'a ProtoCurrent_CSL'

    call PCC % InitializeTemplate_CSL &
           ( C, nValues, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  function ProtoCurrent ( PCC ) result ( PC )

    class ( ProtoCurrent_CSL_Form ), intent ( in ), target :: &
      PCC
    class ( ProtoCurrentForm ), pointer :: &
      PC
      
    class ( VariableGroupForm ), pointer :: &
      Field
    
    Field => PCC % Field
    select type ( Field )
    class is ( ProtoCurrentForm )
    PC => Field
    end select !-- Field

  end function ProtoCurrent


  impure elemental subroutine Finalize ( PCC )

    type ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      PCC

    call PCC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC, NameOption )

    class ( ProtoCurrent_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate ( ProtoCurrentForm :: FC % Field )
    allocate ( FC % FieldOutput )

    select type ( PC => FC % Field )
    class is ( ProtoCurrentForm )
      call PC % Initialize &
             ( FC % VelocityUnit, FC % nValues, NameOption = NameOption )
      call PC % SetOutput ( FC % FieldOutput )
    end select !-- F

  end subroutine SetField


end module ProtoCurrent_CSL__Form
