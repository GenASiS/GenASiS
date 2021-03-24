module FieldSet_MH__Form

  !-- FieldSet_ManifoldHeader__Form

  use Basics
  use Manifold_H__Form

  implicit none
  private

  type, public :: FieldSet_MH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iFieldSet    = 0, &
      nStreams     = 0
    logical ( KDL ) :: &
      Pinned
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
    class ( Manifold_H_Form ), pointer :: &
      Manifold => null ( )
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, public, pass :: &
      Show => Show_FSM
    final :: &
      Finalize
  end type FieldSet_MH_Form

  type, public :: FieldSet_MH_Pointer
    class ( FieldSet_MH_Form ), pointer :: &
      Pointer => null ( )
  end type FieldSet_MH_Pointer


contains


  subroutine Initialize_H &
               ( FSM, M, NameShort, PinnedOption, IgnorabilityOption )

    class ( FieldSet_MH_Form ), intent ( inout ) :: &
      FSM
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      PinnedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FSM % IGNORABILITY  =  M % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FSM % IGNORABILITY  =  IgnorabilityOption

    if ( FSM % Type  ==  '' ) &
      FSM % Type  =  'a FieldSet_M' 
    
    FSM % Pinned  =  .false.
    if ( present ( PinnedOption ) ) &
      FSM % Pinned  =  PinnedOption
    
    FSM % Name  =  trim ( NameShort ) // '_' // trim ( M % Name ) 

    call Show ( 'Initializing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
    FSM % NameShort  =  NameShort

      M % nFieldSets  =  M % nFieldSets  +  1
    FSM % iFieldSet   =  M % nFieldSets

    FSM % Manifold  =>  M

  end subroutine Initialize_H


  subroutine Show_FSM ( FSM )

    class ( FieldSet_MH_Form ), intent ( in ) :: &
      FSM

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FSM % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSM % IGNORABILITY )

    call Show ( FSM % Name,      'Name',      FSM % IGNORABILITY )
    call Show ( FSM % NameShort, 'NameShort', FSM % IGNORABILITY )
    call Show ( FSM % iFieldSet, 'iFieldSet', FSM % IGNORABILITY )
    call Show ( FSM % nStreams,  'nStreams',  FSM % IGNORABILITY )

  end subroutine Show_FSM


  impure elemental subroutine Finalize ( FSM )

    type ( FieldSet_MH_Form ), intent ( inout ) :: &
      FSM

    nullify ( FSM % Manifold )

    call Show ( 'Finalizing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_MH__Form
