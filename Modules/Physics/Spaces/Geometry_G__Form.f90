!-- Geometry_G represents Galilean geometry.

module Geometry_G__Form

  !-- Geometry_Galilean__Form

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_GALILEAN = 0

  type, public, extends ( GeometryFlatForm ) :: Geometry_G_Form
    integer ( KDI ) :: &
      N_FIELDS_GALILEAN = N_FIELDS_GALILEAN
  contains
    procedure, public, pass :: &
      InitializeAllocate_G
    final :: &
      Finalize
  end type Geometry_G_Form


contains


  subroutine InitializeAllocate_G &
               ( G, CoordinateSystem, CoordinateUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 PinnedOption, UnitOption, VectorIndicesOption )

    class ( Geometry_G_Form ), intent ( inout ) :: &
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
      ClearOption, &
      PinnedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    if ( G % Type == '' ) &
      G % Type = 'a Geometry_G'

    call G % GeometryFlatForm % Initialize &
           ( CoordinateSystem, CoordinateUnit, nValues, &
             VariableOption, VectorOption, NameOption, ClearOption, &
             PinnedOption, UnitOption, VectorIndicesOption )

  end subroutine InitializeAllocate_G


  impure elemental subroutine Finalize ( G )

    type ( Geometry_G_Form ), intent ( inout ) :: &
      G

    !-- Empty routine to trigger finalization of parent type

  end subroutine Finalize


end module Geometry_G__Form
