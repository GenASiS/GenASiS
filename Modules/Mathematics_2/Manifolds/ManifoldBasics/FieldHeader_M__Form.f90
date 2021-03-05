!-- FieldHeader_M__Form handles metadata for a set of related fields on 
!   a Manifold.

module FieldHeader_M__Form

  !-- FieldHeader_Manifold__Form

  use Basics
  use ManifoldHeader_Form

  implicit none
  private

  type, public :: FieldHeader_M_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    logical ( KDL ) :: &
      Pinned
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
    class ( ManifoldHeaderForm ), pointer :: &
      Manifold => null ( )
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    final :: &
      Finalize
  end type FieldHeader_M_Form


contains


  subroutine Initialize_H &
               ( FM, M, NameShort, PinnedOption, IgnorabilityOption )

    class ( FieldHeader_M_Form ), intent ( inout ) :: &
      FM
    class ( ManifoldHeaderForm ), intent ( in ), target :: &
      M
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      PinnedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FM % IGNORABILITY = M % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FM % IGNORABILITY = IgnorabilityOption

    if ( FM % Type == '' ) &
      FM % Type = 'a Field_M' 
    
    FM % Pinned = .false.
    if ( present ( PinnedOption ) ) &
      FM % Pinned = PinnedOption
    
    FM % Name = trim ( NameShort ) // '_' // trim ( M % Name ) 

    call Show ( 'Initializing ' // trim ( FM % Type ), FM % IGNORABILITY )
    call Show ( FM % Name, 'Name', FM % IGNORABILITY )
   
    FM % NameShort = NameShort
    call Show ( FM % NameShort, 'NameShort', FM % IGNORABILITY )

    FM % Manifold => M

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FM )

    type ( FieldHeader_M_Form ), intent ( inout ) :: &
      FM

    nullify ( FM % Manifold )

    call Show ( 'Finalizing ' // trim ( FM % Type ), FM % IGNORABILITY )
    call Show ( FM % Name, 'Name', FM % IGNORABILITY )
   
  end subroutine Finalize


end module FieldHeader_M__Form
