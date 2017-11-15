!-- Geometry_CSL represents the geometry on a single-level Chart.

module Geometry_CSL__Form

  !-- GeometryFlat_ChartSingleLevel_Form

  use Basics
  use Mathematics
  use Geometry_G__Form

  implicit none
  private

  type, public, extends ( GeometryFlat_CSL_Form ) :: Geometry_CSL_Form
  contains
    procedure, public, pass :: &
      InitializeType
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Geometry_CSL_Form

contains


  subroutine InitializeType &
               ( GC, C, NameShort, GeometryType, nValues, IgnorabilityOption )

    class ( Geometry_CSL_Form ), intent ( inout ) :: &
      GC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      GeometryType
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( GC % Type == '' ) &
      GC % Type = 'a Geometry_CSL'

    GC % GeometryType = GeometryType    

    call GC % Initialize ( C, NameShort, nValues, IgnorabilityOption )

  end subroutine InitializeType


  impure elemental subroutine Finalize ( GC )

    type ( Geometry_CSL_Form ), intent ( inout ) :: &
      GC

    !-- Trigger finalization of parent

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( Geometry_CSL_Form ), intent ( inout ) :: &
      FC

    select case ( trim ( FC % GeometryType ) )
    case ( 'GALILEAN' )
      allocate ( GeometryFlatForm :: FC % Field )
    case default
      call Show ( 'GeometryType not recognized', CONSOLE % ERROR )
      call Show ( FC % GeometryType, 'GeometryType', CONSOLE % ERROR )
      call Show ( 'Geometry_CSL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

    allocate ( FC % FieldOutput )

    select type ( G => FC % Field )
    class is ( Geometry_G_Form )
      call G % Initialize &
             ( FC % Chart % CoordinateSystem, FC % Chart % CoordinateUnit, &
               FC % nValues, NameOption = FC % NameShort )
      call G % SetOutput ( FC % FieldOutput )
    end select !-- F

  end subroutine SetField


end module Geometry_CSL__Form