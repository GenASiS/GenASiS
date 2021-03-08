module FieldSet_MH__Form

  !-- FieldSet_ManifoldHeader__Form

  use Basics
  use ManifoldHeader_Form

  implicit none
  private

  type, public :: FieldSet_MH_Form
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
  end type FieldSet_MH_Form

  type, public :: FieldSet_MH_Pointer
    class ( FieldSet_MH_Form ), pointer :: &
      Pointer => null ( )
  end type FieldSet_MH_Pointer


contains


  subroutine Initialize_H &
               ( FM, M, NameShort, PinnedOption, IgnorabilityOption )

    class ( FieldSet_MH_Form ), intent ( inout ) :: &
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
      FM % Type = 'a FieldSet_M' 
    
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

    type ( FieldSet_MH_Form ), intent ( inout ) :: &
      FM

    nullify ( FM % Manifold )

    call Show ( 'Finalizing ' // trim ( FM % Type ), FM % IGNORABILITY )
    call Show ( FM % Name, 'Name', FM % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_MH__Form
