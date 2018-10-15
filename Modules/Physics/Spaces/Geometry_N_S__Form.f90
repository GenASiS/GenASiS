!-- Geometry_N_S represents Newtonian geometry, plus an overloading of
!   ComputeReconstruction to add reconstruction of gradient of Phi for use
!   in the divergence of a gravitational stress as an alternative to a 
!   gravitational force source term.

module Geometry_N_S__Form

  !-- Geometry_Newtonian_Stress__Form

  use Basics
  use Mathematics
  use Geometry_N__Form

  implicit none
  private

  type, public, extends ( Geometry_N_Form ) :: Geometry_N_S_Form
  contains
    procedure, public, pass :: &
      InitializeAllocate_G
    final :: &
      Finalize
  end type Geometry_N_S_Form

contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( Geometry_N_S_Form ), intent ( inout ) :: &
      G
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      CoordinateUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    if ( G % Type == '' ) &
      G % Type = 'a Geometry_N_S'

    call G % Geometry_N_Form % Initialize &
           ( CoordinateSystem, CoordinateUnit, nValues, &
             VariableOption, VectorOption, NameOption, ClearOption, &
             UnitOption, VectorIndicesOption )

  end subroutine InitializeAllocate_G


  impure elemental subroutine Finalize ( G )

    type ( Geometry_N_S_Form ), intent ( inout ) :: &
      G

    !-- Empty routine to trigger finalization of parent type

  end subroutine Finalize


end module Geometry_N_S__Form
