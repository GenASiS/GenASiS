module FieldSet_CH__Form

  !-- FieldSet_ChartHeader__Form

  use Basics
  use ManifoldBasics
  use Chart_H__Form

  implicit none
  private

  type, public :: FieldSet_CH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
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
               ( FC, FSM, C, NameShort, PinnedOption, IgnorabilityOption )

    class ( FieldSet_CH_Form ), intent ( inout ) :: &
      FC
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

    FC % IGNORABILITY = C % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FC % IGNORABILITY = IgnorabilityOption

    if ( FC % Type == '' ) &
      FC % Type = 'a FieldSet_C' 
    
    FC % Pinned = .false.
    if ( present ( PinnedOption ) ) &
      FC % Pinned = PinnedOption
    
    FC % Name = trim ( NameShort ) // '_' // trim ( C % Name ) 

    call Show ( 'Initializing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
    FC % NameShort = NameShort
    call Show ( FC % NameShort, 'NameShort', FC % IGNORABILITY )

    FC % Chart       =>    C
    FC % FieldSet_M  =>  FSM 

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FC )

    type ( FieldSet_CH_Form ), intent ( inout ) :: &
      FC

    nullify ( FC % FieldSet_M )
    nullify ( FC % Chart )

    call Show ( 'Finalizing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_CH__Form
