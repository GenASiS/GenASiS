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
      DeviceMemory, &
      PinnedMemory, &
      DevicesCommunicate
    character ( LDL ) :: &
      Name = '', &
      Type = ''
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
               ( FSM, M, Name, DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption )

    class ( FieldSet_MH_Form ), intent ( inout ) :: &
      FSM
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    character ( * ), intent ( in ) :: &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption

    FSM % IGNORABILITY  =  M % IGNORABILITY

    if ( FSM % Type  ==  '' ) &
      FSM % Type  =  'a FieldSet_M' 
    
    FSM % DeviceMemory  =  .false.
    if ( present ( DeviceMemoryOption ) ) &
      FSM % DeviceMemory  =  DeviceMemoryOption
    
    FSM % PinnedMemory  =  .false.
    if ( present ( PinnedMemoryOption ) ) &
      FSM % PinnedMemory  =  PinnedMemoryOption
    
    FSM % DevicesCommunicate  =  .false.
    if ( present ( DevicesCommunicateOption ) )  &
      FSM % DevicesCommunicate  =  DevicesCommunicateOption  

    FSM % Name  =  Name

    call Show ( 'Initializing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
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

    associate ( M  =>  FSM % Manifold )
    call Show ( FSM % Name,      'Name',      FSM % IGNORABILITY )
    call Show (   M % Name,      'Manifold',  FSM % IGNORABILITY )
    call Show ( FSM % iFieldSet, 'iFieldSet', FSM % IGNORABILITY )
    end associate  !-- M

    call Show ( FSM % DeviceMemory,       'DeviceMemory', &
                FSM % IGNORABILITY )
    call Show ( FSM % PinnedMemory,       'PinnedMemory', &
                FSM % IGNORABILITY )
    call Show ( FSM % DevicesCommunicate, 'DevicesCommunicate', &
                FSM % IGNORABILITY )

  end subroutine Show_FSM


  impure elemental subroutine Finalize ( FSM )

    type ( FieldSet_MH_Form ), intent ( inout ) :: &
      FSM

    nullify ( FSM % Manifold )

    call Show ( 'Finalizing ' // trim ( FSM % Type ), FSM % IGNORABILITY )
    call Show ( FSM % Name, 'Name', FSM % IGNORABILITY )
   
  end subroutine Finalize


end module FieldSet_MH__Form
