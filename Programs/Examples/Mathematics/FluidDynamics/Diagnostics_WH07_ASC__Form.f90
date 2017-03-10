module Diagnostics_WH07_ASC__Form

  use Basics
  use Mathematics
  use Diagnostics_WH07__Form
  use Diagnostics_WH07_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Diagnostics_WH07_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Diagnostics_WH07_CSL
    generic, public :: &
      Diagnostics_WH07 => Diagnostics_WH07_CSL
    procedure, private, pass :: &
      SetField
  end type Diagnostics_WH07_ASC_Form

contains


  subroutine Initialize ( DA, A, NameShort )

    class ( Diagnostics_WH07_ASC_Form ), intent ( inout ) :: &
      DA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      NameShort

    call DA % InitializeTemplate_ASC ( A, NameShort )

  end subroutine Initialize


  function Diagnostics_WH07_CSL ( DA ) result ( D )

    class ( Diagnostics_WH07_ASC_Form ), intent ( in ) :: &
      DA
    class ( Diagnostics_WH07_Form ), pointer :: &
      D

    select type ( DC => DA % Chart )
    class is ( Diagnostics_WH07_CSL_Form )
      D => DC % Diagnostics ( )
    class default
      call Show ( 'Diagnostics Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Diagnostics_WH07_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Diagnostics_WH07_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Diagnostics_WH07_CSL


  subroutine SetField ( FA )

    class ( Diagnostics_WH07_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Diagnostics_WH07_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Diagnostics_WH07_CSL_Form )
      call FC % Initialize ( C, FA % NameShort, nValues )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Diagnostics_WH07_ASC__Form
