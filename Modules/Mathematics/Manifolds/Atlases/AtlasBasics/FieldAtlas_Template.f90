!-- FieldAtlas is a template for a set of fields on a Atlas.

module FieldAtlas_Template

  use Basics
  use AtlasHeader_Form

  implicit none
  private

  type, public, abstract :: FieldAtlasTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameOutput = ''
    class ( AtlasHeaderForm ), pointer :: &
      Atlas => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type FieldAtlasTemplate

contains


  subroutine InitializeTemplate ( FA, A, NameOutputOption )

    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      FA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    FA % IGNORABILITY = A % IGNORABILITY

    if ( FA % Type == '' ) &
      FA % Type = 'a FieldAtlas' 

    call Split ( FA % Type, ' ', TypeWord )
    FA % Name = trim ( TypeWord ( 2 ) ) // '_' // trim ( A % Name ) 

    call Show ( 'Initializing ' // trim ( FA % Type ), FA % IGNORABILITY )
    call Show ( FA % Name, 'Name', FA % IGNORABILITY )
   
    if ( present ( NameOutputOption ) ) then
      FA % NameOutput = NameOutputOption
      call Show ( FA % NameOutput, 'NameOutput', FA % IGNORABILITY )
    end if

    FA % Atlas => A

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( FA )

    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      FA

    nullify ( FA % Atlas )

    call Show ( 'Finalizing ' // trim ( FA % Type ), FA % IGNORABILITY )
    call Show ( FA % Name, 'Name', FA % IGNORABILITY )
   
  end subroutine FinalizeTemplate


end module FieldAtlas_Template
