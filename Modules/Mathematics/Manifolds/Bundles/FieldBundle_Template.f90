!-- FieldBundle is a template for a set of fields on a Bundle.

module FieldBundle_Template

  use Basics
  use BundleHeader_Form

  implicit none
  private

  type, public, abstract :: FieldBundleTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0
    logical ( KDL ) :: &
      UsePinnedMemory_S, &
      UsePinnedMemory_F
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
    class ( BundleHeaderForm ), pointer :: &
      Bundle => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type FieldBundleTemplate

    abstract interface 
      subroutine SF ( FB )
        import FieldBundleTemplate
        class ( FieldBundleTemplate ), intent ( inout ) :: &
          FB
      end subroutine
    end interface

  ! type, public :: FieldBundleElementForm
  !   class ( FieldBundleTemplate ), allocatable :: &
  !     Element
  ! contains
  !   final :: &
  !     FinalizeElement
  ! end type FieldBundleElementForm

contains


  subroutine InitializeTemplate &
               ( FB, B, NameShort, UsePinnedMemory_S_Option, &
                 UsePinnedMemory_F_Option, IgnorabilityOption )

    class ( FieldBundleTemplate ), intent ( inout ) :: &
      FB
    class ( BundleHeaderForm ), intent ( in ), target :: &
      B
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemory_S_Option, &
      UsePinnedMemory_F_Option
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FB % IGNORABILITY = B % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FB % IGNORABILITY = IgnorabilityOption

    if ( FB % Type == '' ) &
      FB % Type = 'a FieldBundle' 
      
    FB % UsePinnedMemory_S = .false.
    FB % UsePinnedMemory_F = .false.
    if ( present ( UsePinnedMemory_S_Option ) ) &
      FB % UsePinnedMemory_S = UsePinnedMemory_S_Option
    if ( present ( UsePinnedMemory_F_Option ) ) &
      FB % UsePinnedMemory_F = UsePinnedMemory_F_Option

    FB % Name = trim ( NameShort ) // '_' // B % Name 

    call Show ( 'Initializing ' // trim ( FB % Type ), FB % IGNORABILITY )
    call Show ( FB % Name, 'Name', FB % IGNORABILITY )
   
    FB % NameShort = NameShort
    call Show ( FB % NameShort, 'NameShort', FB % IGNORABILITY )

    FB % Bundle => B

    call FB % SetField ( )

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
