module Diagnostics_WH07_CSL__Form

  use Basics
  use Mathematics
  use Diagnostics_WH07__Form
  
  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Diagnostics_WH07_CSL_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Diagnostics
    procedure, private, pass :: &
      SetField
  end type Diagnostics_WH07_CSL_Form

contains


  subroutine Initialize ( DC, C, NameShort, nValues )

    class ( Diagnostics_WH07_CSL_Form ), intent ( inout ) :: &
      DC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nValues

    call DC % InitializeTemplate_CSL ( C, NameShort, nValues )

  end subroutine Initialize


  function Diagnostics ( DC ) result ( D )

    class ( Diagnostics_WH07_CSL_Form ), intent ( in ), target :: &
      DC
    class ( Diagnostics_WH07_Form ), pointer :: &
      D
      
    class ( VariableGroupForm ), pointer :: &
      Field

    D => null ( )

    Field => DC % Field
    select type ( Field )
    class is ( Diagnostics_WH07_Form )
    D => Field
    end select !-- Field

  end function Diagnostics


  subroutine SetField ( FC )

    class ( Diagnostics_WH07_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % FieldOutput )
    allocate ( Diagnostics_WH07_Form :: FC % Field )
    select type ( D => FC % Field )
    class is ( Diagnostics_WH07_Form )
      call D % Initialize ( FC % NameShort, FC % nValues )
      call D % SetOutput ( FC % FieldOutput )
    end select !-- D

  end subroutine SetField

end module Diagnostics_WH07_CSL__Form
