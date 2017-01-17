!-- FieldBundle is a template for a set of fields on a Bundle.

module FieldBundle_Template

  use Basics
  use BundleHeader_Form

  implicit none
  private

  type, public, abstract :: FieldBundleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameOutput = ''
    class ( BundleHeaderForm ), pointer :: &
      Bundle => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type FieldBundleTemplate

  ! type, public :: FieldBundleElementForm
  !   class ( FieldBundleTemplate ), allocatable :: &
  !     Element
  ! contains
  !   final :: &
  !     FinalizeElement
  ! end type FieldBundleElementForm

contains


  subroutine InitializeTemplate ( FB, B, NameOutputOption )

    class ( FieldBundleTemplate ), intent ( inout ) :: &
      FB
    class ( BundleHeaderForm ), intent ( in ), target :: &
      B
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    FB % IGNORABILITY = B % IGNORABILITY

    if ( FB % Type == '' ) &
      FB % Type = 'a FieldBundle' 

    call Split ( FB % Type, ' ', TypeWord )
    FB % Name = trim ( B % Name ) // '_' // trim ( TypeWord ( 2 ) ) 

    call Show ( 'Initializing ' // trim ( FB % Type ), FB % IGNORABILITY )
    call Show ( FB % Name, 'Name', FB % IGNORABILITY )
   
    if ( present ( NameOutputOption ) ) then
      FB % NameOutput = NameOutputOption
      call Show ( FB % NameOutput, 'NameOutput', FB % IGNORABILITY )
    end if

    FB % Bundle => B

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( FB )

    class ( FieldBundleTemplate ), intent ( inout ) :: &
      FB

    nullify ( FB % Bundle )

    call Show ( 'Finalizing ' // trim ( FB % Type ), FB % IGNORABILITY )
    call Show ( FB % Name, 'Name', FB % IGNORABILITY )
   
  end subroutine FinalizeTemplate


  ! impure elemental subroutine FinalizeElement ( FBE )
    
  !   type ( FieldBundleElementForm ), intent ( inout ) :: &
  !     FBE

  !   if ( allocated ( FBE % Element ) ) deallocate ( FBE % Element )

  ! end subroutine FinalizeElement


end module FieldBundle_Template
