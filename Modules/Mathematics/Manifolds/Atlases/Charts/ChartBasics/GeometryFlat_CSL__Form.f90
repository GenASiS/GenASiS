!-- GeometryFlat_CSL represents the geometry on a single-level Chart.

module GeometryFlat_CSL__Form

  !-- GeometryFlat_ChartSingleLevel_Form

  use Basics
  use ChartHeader_SL__Form
  use GeometryFlat_Form
  use Field_CSL__Template

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: GeometryFlat_CSL_Form
    character ( LDL ) :: &
      GeometryType = ''
  contains
    procedure, public, pass :: &
      InitializeFlat
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type GeometryFlat_CSL_Form

contains


  subroutine InitializeFlat &
               ( GC, C, NameShort, UsePinnedMemory, nValues, &
                 IgnorabilityOption )

    class ( GeometryFlat_CSL_Form ), intent ( inout ) :: &
      GC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ) :: &
      UsePinnedMemory
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( GC % Type == '' ) &
      GC % Type = 'a GeometryFlat_CSL'

    if ( GC % GeometryType == '' ) &
      GC % GeometryType = 'FLAT'    

    call GC % InitializeTemplate_CSL &
           ( C, NameShort, UsePinnedMemory, nValues, IgnorabilityOption )

  end subroutine InitializeFlat


  impure elemental subroutine Finalize ( GC )

    type ( GeometryFlat_CSL_Form ), intent ( inout ) :: &
      GC

    call GC % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( GeometryFlat_CSL_Form ), intent ( inout ) :: &
      FC
                      
    if ( .not. allocated ( FC % Field ) ) then
      select case ( trim ( FC % GeometryType ) )
      case ( 'FLAT' )
        allocate ( GeometryFlatForm :: FC % Field )
      case default
        call Show ( 'GeometryType not recognized', CONSOLE % ERROR )
        call Show ( FC % GeometryType, 'GeometryType', CONSOLE % ERROR )
        call Show ( 'GeometryFlat_CSL__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- FluidType
    end if !-- alloc Field % Element

    allocate ( FC % FieldOutput )

    select type ( G => FC % Field )
    class is ( GeometryFlatForm )
      call G % Initialize &
             ( FC % Chart % CoordinateSystem, FC % Chart % CoordinateUnit, &
               FC % nValues, NameOption = FC % NameShort )
      call G % SetOutput ( FC % FieldOutput )
    end select !-- F

  end subroutine SetField


end module GeometryFlat_CSL__Form
