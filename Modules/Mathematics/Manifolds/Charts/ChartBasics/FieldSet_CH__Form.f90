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
      iFieldSet    = 0, &
      nStreams     = 0
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
    procedure, public, pass :: &
      Show => Show_FSC
    final :: &
      Finalize
  end type FieldSet_CH_Form

  type, public :: FieldSet_CH_Pointer
    class ( FieldSet_CH_Form ), pointer :: &
      Pointer => null ( )
  end type FieldSet_CH_Pointer


contains


  subroutine Initialize_H ( FSC, C, FSM )

    class ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( inout ), target :: &
      C
    class ( FieldSet_MH_Form ), intent ( in ), target :: &
      FSM

    FSC % IGNORABILITY  =  C % IGNORABILITY

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a FieldSet_C' 
    
    FSC % Pinned  =  FSM % Pinned
    
    FSC % Name  =  trim ( FSM % NameShort ) // '_' // trim ( C % Name ) 

    call Show ( 'Initializing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
      C % nFieldSets  =    C % Manifold % nFieldSets
    FSC % iFieldSet   =  FSM % iFieldSet

    FSC % NameShort  =  FSM % NameShort

    FSC % Chart       =>    C
    FSC % FieldSet_M  =>  FSM 

  end subroutine Initialize_H


  subroutine Show_FSC ( FSC )

    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( FSC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSC % IGNORABILITY )

    call Show ( FSC % Name,      'Name',      FSC % IGNORABILITY )
    call Show ( FSC % NameShort, 'NameShort', FSC % IGNORABILITY )
    call Show ( FSC % iFieldSet, 'iFieldSet', FSC % IGNORABILITY )

  end subroutine Show_FSC


  impure elemental subroutine Finalize ( FSC )

    type ( FieldSet_CH_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % FieldSet_M )
    nullify ( FSC % Chart )

    call Show ( 'Finalizing ' // trim ( FSC % Type ), FSC % IGNORABILITY )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_CH__Form
