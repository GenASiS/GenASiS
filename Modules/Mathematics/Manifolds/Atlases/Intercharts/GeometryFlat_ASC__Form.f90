!-- GeometryFlat_ASC represents the geometry on a single-chart Atlas.

module GeometryFlat_ASC__Form

  !-- GeometryFlat_AtlasSingleChart_Form

  use Basics
  use AtlasBasics
  use Charts
  use Field_ASC__Template
  use Atlas_SC__Template

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: GeometryFlat_ASC_Form
    character ( LDL ) :: &
      GeometryType = ''
  contains
    procedure, private, pass :: &
      InitializeFlat
    generic, public :: &
      Initialize => InitializeFlat
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type GeometryFlat_ASC_Form

contains


  subroutine InitializeFlat ( GA, A, NameShortOption, IgnorabilityOption )

    class ( GeometryFlat_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      NameShort

    if ( GA % Type == '' ) &
      GA % Type = 'a GeometryFlat_ASC'

    if ( GA % GeometryType == '' ) &
      GA % GeometryType = 'FLAT'    

    NameShort = 'Geometry'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call GA % InitializeTemplate_ASC ( A, NameShort, IgnorabilityOption )

    call Show ( GA % GeometryType, 'GeometryType', GA % IGNORABILITY )

  end subroutine InitializeFlat


  impure elemental subroutine Finalize ( GA )

    type ( GeometryFlat_ASC_Form ), intent ( inout ) :: &
      GA

    call GA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( GeometryFlat_ASC_Form ), intent ( inout ) :: &
      FA

    if ( .not. allocated ( FA % Chart ) ) then
      select case ( trim ( FA % GeometryType ) )
      case ( 'FLAT' )

        select type ( A => FA % Atlas )
        class is ( Atlas_SC_Template )

        select type ( C => A % Chart )
        class is ( Chart_SL_Template )

          allocate ( GeometryFlat_CSL_Form :: FA % Chart )

          select type ( GC => FA % Chart )
          class is ( GeometryFlat_CSL_Form )
            associate ( nValues => C % nProperCells + C % nGhostCells )
            call GC % Initialize &
                   ( C, FA % NameShort, nValues, &
                     IgnorabilityOption = FA % IGNORABILITY )
            end associate !-- nValues
          end select !-- GC

        end select !-- C
        end select !-- A

      case default
        call Show ( 'GeometryType not recognized', CONSOLE % ERROR )
        call Show ( FA % GeometryType, 'GeometryType', CONSOLE % ERROR )
        call Show ( 'GeometryFlat_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetField', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- GeometryType
    end if !-- alloc Field % Element

  end subroutine SetField


end module GeometryFlat_ASC__Form
