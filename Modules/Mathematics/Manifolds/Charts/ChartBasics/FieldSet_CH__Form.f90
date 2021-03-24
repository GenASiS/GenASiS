module FieldSet_CH__Form

  !-- FieldSet_ChartHeader__Form

  use Basics
  use ManifoldBasics
  use Chart_H__Form

  implicit none
  private

  type, public :: FieldSet_CH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iFieldSet    = 0
    logical ( KDL ) :: &
      Pinned = .false.
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
    class ( FieldSet_MH_Form ), pointer :: &
      FieldSet_M => null ( )
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    final :: &
      Finalize
  end type FieldSet_CH_Form

  type, public :: FieldSet_CH_Pointer
    class ( FieldSet_CH_Form ), pointer :: &
      Pointer => null ( )
  end type FieldSet_CH_Pointer


contains


  subroutine Initialize_H &
               ( FSC, FSM, C, NameShort, PinnedOption, IgnorabilityOption )

    class ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC
    class ( FieldSet_MH_Form ), intent ( in ), target :: &
      FSM
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      PinnedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FSC % IGNORABILITY  =  C % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FSC % IGNORABILITY  =  IgnorabilityOption

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a FieldSet_C' 
    
    FSC % Pinned  =  .false.
    if ( present ( PinnedOption ) ) &
      FSC % Pinned  =  PinnedOption
    
    FSC % Name  =  trim ( NameShort ) // '_' // trim ( C % Name ) 

    call Show ( 'Initializing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
    FSC % NameShort  =  NameShort
    call Show ( FSC % NameShort, 'NameShort', FSC % IGNORABILITY )

      C % nFieldSets  =    M % nFieldSets
    FSC % iFieldSet   =  FSM % iFieldSet
    call Show ( FSC % iFieldSet, 'iFieldSet' )

    FSC % Chart       =>    C
    FSC % FieldSet_M  =>  FSM 

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FSC )

    type ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % FieldSet_M )
    nullify ( FSC % Chart )

    call Show ( 'Finalizing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_CH__Form
