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
      NameShort = ''
    class ( AtlasHeaderForm ), pointer :: &
      Atlas => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type FieldAtlasTemplate

    abstract interface 
      subroutine SF ( FA )
        import FieldAtlasTemplate
        class ( FieldAtlasTemplate ), intent ( inout ) :: &
          FA
      end subroutine
    end interface

  type, public :: FieldAtlasPointer
    class ( FieldAtlasTemplate ), pointer :: &
      Pointer => null ( )
  end type FieldAtlasPointer

  type, public :: FieldAtlasElement
    class ( FieldAtlasTemplate ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type FieldAtlasElement

  type, public :: FieldAtlas_1D_Form
    type ( FieldAtlasElement ), dimension ( : ), allocatable :: &
      Atlas
  contains
    procedure, public, pass :: &
      Initialize => Initialize_1D
    final :: &
      Finalize_1D
  end type FieldAtlas_1D_Form

contains


  subroutine InitializeTemplate ( FA, A, NameShort, IgnorabilityOption )

    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      FA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FA % IGNORABILITY = A % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FA % IGNORABILITY = IgnorabilityOption

    if ( FA % Type == '' ) &
      FA % Type = 'a FieldAtlas' 

    FA % Name = trim ( NameShort ) // '_' // trim ( A % Name ) 

    call Show ( 'Initializing ' // trim ( FA % Type ), FA % IGNORABILITY )
    call Show ( FA % Name, 'Name', FA % IGNORABILITY )
   
    FA % NameShort = NameShort
    call Show ( FA % NameShort, 'NameShort', FA % IGNORABILITY )

    FA % Atlas => A

    call FA % SetField ( )

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( FA )

    class ( FieldAtlasTemplate ), intent ( inout ) :: &
      FA

    nullify ( FA % Atlas )

    call Show ( 'Finalizing ' // trim ( FA % Type ), FA % IGNORABILITY )
    call Show ( FA % Name, 'Name', FA % IGNORABILITY )
   
  end subroutine FinalizeTemplate


  impure elemental subroutine FinalizeElement ( FAE )
    
    type ( FieldAtlasElement ), intent ( inout ) :: &
      FAE

    if ( allocated ( FAE % Element ) ) &
      deallocate ( FAE % Element )

  end subroutine FinalizeElement


  subroutine Initialize_1D ( FA_1D, nElements )

    class ( FieldAtlas_1D_Form ), intent ( inout ) :: &
      FA_1D
    integer ( KDI ), intent ( in ) :: &
      nElements

    allocate ( FA_1D % Atlas ( nElements ) )

  end subroutine Initialize_1D


  impure elemental subroutine Finalize_1D ( FA_1D )

    type ( FieldAtlas_1D_Form ), intent ( inout ) :: &
      FA_1D
    
    if ( allocated ( FA_1D % Atlas ) ) &
      deallocate ( FA_1D % Atlas )

  end subroutine Finalize_1D


end module FieldAtlas_Template
