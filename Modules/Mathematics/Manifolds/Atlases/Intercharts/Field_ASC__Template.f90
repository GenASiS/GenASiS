!-- Field_ASC is a template for a set of fields on an Atlas_SC.

module Field_ASC__Template

  !-- Field_AtlasSingleChart_Template

  use Basics
  use AtlasBasics
  use Charts

  implicit none
  private

  type, public, extends ( FieldAtlasTemplate ), abstract :: Field_ASC_Template
    class ( FieldChartTemplate ), allocatable :: &
      Chart
  contains
    procedure, public, pass :: &
      InitializeTemplate_ASC
    procedure, public, pass :: &
      FinalizeTemplate_ASC
  end type Field_ASC_Template

contains


  subroutine InitializeTemplate_ASC ( FA, A, NameShort, IgnorabilityOption )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FA % Type == '' ) &
      FA % Type = 'a Field_ASC' 

    call FA % InitializeTemplate ( A, NameShort, IgnorabilityOption )

  end subroutine InitializeTemplate_ASC


  impure elemental subroutine FinalizeTemplate_ASC ( FA )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA

    if ( allocated ( FA % Chart ) ) &
      deallocate ( FA % Chart )

    call FA % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_ASC


end module Field_ASC__Template
