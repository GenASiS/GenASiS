!-- Geometry_ASC represents the geometry on a position space represented
!   by a single chart.

module Geometry_ASC__Form

  !-- Geometry_AtlasSingleChart_Form

  use Basics
  use Mathematics
  use Geometry_CSL__Form

  implicit none
  private

  type, public, extends ( GeometryFlat_ASC_Form ) :: Geometry_ASC_Form
  contains
    procedure, public, pass :: &
      InitializeType
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Geometry_ASC_Form

contains


  subroutine InitializeType &
               ( GA, A, GeometryType, NameShortOption, IgnorabilityOption )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      GeometryType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      NameShort

    if ( GA % Type == '' ) &
      GA % Type = 'a Geometry_ASC'

    GA % GeometryType = GeometryType    

    call GA % Initialize ( A, NameShortOption, IgnorabilityOption )

  end subroutine InitializeType


  impure elemental subroutine Finalize ( GA )

    type ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA

    !-- Trigger finalization of parent

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Geometry_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Geometry_CSL_Form )
      call FC % InitializeType &
             ( C, FA % NameShort, FA % GeometryType, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Geometry_ASC__Form
