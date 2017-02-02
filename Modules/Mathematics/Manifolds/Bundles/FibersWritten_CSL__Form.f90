!-- FibersWritten_CSL flags the Fibers written on a single-level Chart.

module FibersWritten_CSL__Form

  !-- FibersWritten_ChartSingleLevel

  use Basics
  use Atlases

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: FibersWritten_CSL_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type FibersWritten_CSL_Form

contains


  subroutine Initialize ( FWC, C, nValues, NameOutputOption )

    class ( FibersWritten_CSL_Form ), intent ( inout ) :: &
      FWC
    class ( Chart_SL_Template ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
 
    if ( FWC % Type == '' ) &
      FWC % Type = 'a FibersWritten_CSL'

    call FWC % InitializeTemplate_CSL &
           ( C, nValues, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FWC )

    type ( FibersWritten_CSL_Form ), intent ( inout ) :: &
      FWC

    call FWC % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetField ( FC, NameOption )

    class ( FibersWritten_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate &
      ( FC % Field, &
        FC % FieldOutput )
    associate &
      ( FW  => FC % Field, &
        FWO => FC % FieldOutput )

    call FW % Initialize &
           ( [ FC % nValues, 1 ], VariableOption = [ 'Flag' ], &
             NameOption = 'FibersWritten', ClearOption = .true. )
    call FWO % Initialize ( FW )

    end associate !-- FW, FWO

  end subroutine SetField


end module FibersWritten_CSL__Form
