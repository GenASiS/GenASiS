module WoosleyHeger_07_A__Form

  !-- WoosleyHeger_07_Adiabatic__Form

  use GenASiS
  use WoosleyHeger_07__Template

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Template ) :: WoosleyHeger_07_A_Form
  contains
    procedure, public, pass :: &
      Initialize
  end type WoosleyHeger_07_A_Form

contains


  subroutine Initialize ( WH, Name )

    class ( WoosleyHeger_07_A_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

  end subroutine Initialize


end module WoosleyHeger_07_A__Form
