module Interactions_B85__Form

  !-- Interactions_Bruenn_1985__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_B85_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      ComputeTimeScale
  end type Interactions_B85_Form

contains


  subroutine InitializeAllocate_I &
               ( I, MomentsType, Units, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_B85_Form ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      MomentsType
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    !-- FIXME: fill in

  end subroutine InitializeAllocate_I


  subroutine Compute ( I, R )

    class ( Interactions_B85_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    !-- FIXME: fill in

  end subroutine Compute


  subroutine ComputeTimeScale ( I, R )

    class ( Interactions_B85_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    !-- FIXME: fill in

  end subroutine ComputeTimeScale


end module Interactions_B85__Form
