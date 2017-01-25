module Interactions_F__Form

  !-- Compute does nothing, leaving interactions as set at initialization.
  
  use Basics
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_F_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      Compute
    generic, public :: &
      Initialize => InitializeAllocate_F
    final :: &
      Finalize
  end type Interactions_F_Form

contains


  subroutine InitializeAllocate_F &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_F_Form ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_F'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, nValues, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_F


  subroutine Compute ( I )

    class ( Interactions_F_Form ), intent ( inout ) :: &
      I

    !-- Compute does nothing, leaving interactions as set at initialization.
  
  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_F_Form ), intent ( inout ) :: &
      I

    call I % FinalizeTemplate ( )

  end subroutine Finalize


end module Interactions_F__Form
