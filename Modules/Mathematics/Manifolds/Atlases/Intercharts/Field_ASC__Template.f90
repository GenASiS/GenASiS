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
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type Field_ASC_Template

    abstract interface 
      subroutine SF ( FA, NameOutputOption )
        import Field_ASC_Template
        class ( Field_ASC_Template ), intent ( inout ) :: &
          FA
        character ( * ), intent ( in ), optional :: &
          NameOutputOption
      end subroutine
    end interface

  type, public :: Field_ASC_Pointer
    class ( Field_ASC_Template ), pointer :: &
      Pointer => null ( )
  end type Field_ASC_Pointer

  type, public :: Field_ASC_ElementForm
    class ( Field_ASC_Template ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type Field_ASC_ElementForm

  type, public :: Field_ASC_1D_Form
    class ( Field_ASC_ElementForm ), dimension ( : ), allocatable :: &
      Atlas
  contains
    procedure, public, pass :: &
      Initialize => Initialize_1D
    final :: &
      Finalize_1D
  end type Field_ASC_1D_Form

contains


  subroutine InitializeTemplate_ASC ( FA, A, NameOutputOption )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    call FA % InitializeTemplate ( A, NameOutputOption = NameOutputOption )

    if ( FA % NameOutput == '' ) then
      call FA % SetField ( )
    else
      call FA % SetField ( NameOutputOption = FA % NameOutput )
    end if

  end subroutine InitializeTemplate_ASC


  impure elemental subroutine FinalizeTemplate_ASC ( FA )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA

    if ( allocated ( FA % Chart ) ) &
      deallocate ( FA % Chart )

    call FA % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_ASC


  impure elemental subroutine FinalizeElement ( FAE )
    
    type ( Field_ASC_ElementForm ), intent ( inout ) :: &
      FAE

    if ( allocated ( FAE % Element ) ) deallocate ( FAE % Element )

  end subroutine FinalizeElement


  subroutine Initialize_1D ( FA_1D, nElements )

    class ( Field_ASC_1D_Form ), intent ( inout ) :: &
      FA_1D
    integer ( KDI ), intent ( in ) :: &
      nElements

    allocate ( FA_1D % Atlas ( nElements ) )

  end subroutine Initialize_1D


  impure elemental subroutine Finalize_1D ( FA_1D )

    type ( Field_ASC_1D_Form ), intent ( inout ) :: &
      FA_1D
    
    if ( allocated ( FA_1D % Atlas ) ) &
      deallocate ( FA_1D % Atlas )

  end subroutine Finalize_1D


end module Field_ASC__Template
