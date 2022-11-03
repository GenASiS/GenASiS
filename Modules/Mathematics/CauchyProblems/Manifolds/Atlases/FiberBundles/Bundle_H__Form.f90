module Bundle_H__Form

  !-- Bundle_Header__Form

  use Basics
  use BaseManifolds

  implicit none
  private

  type, public :: Bundle_H_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( Atlas_H_Form ), allocatable :: &
      Fiber
    class ( Atlas_H_Form ), pointer :: &
      Base => null ( )
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, private, pass :: &
      Show_B
    generic, public :: &
      Show => Show_B
    final :: &
      Finalize
  end type Bundle_H_Form


contains


  subroutine Initialize_H ( B, Base, NameOption )

    class ( Bundle_H_Form ), intent ( inout ) :: &
      B
    class ( Atlas_H_Form ), intent ( in ), target :: &
      Base
    character ( * ), intent ( in ), optional :: &
      NameOption

    B % IGNORABILITY  =  CONSOLE % INFO_1

    if ( B % Type  ==  '' ) &
      B % Type  =  'a Bundle'

    B % Name  =  'Bundle'
    if ( present ( NameOption ) ) &
      B % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( B % Type ), B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )

  end subroutine Initialize_H


  subroutine Show_B ( B )

    class ( Bundle_H_Form ), intent ( in ) :: &
      B

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( B % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )

    call Show ( B % Base % Name, 'Base', B % IGNORABILITY )

    if ( allocated ( B % Fiber ) ) &
      call B % Fiber % Show ( )

  end subroutine Show_B


  impure elemental subroutine Finalize ( B )

    type ( Bundle_H_Form ), intent ( inout ) :: &
      B

    nullify ( B % Base )

    if ( allocated ( B % Fiber ) ) &
      deallocate ( B % Fiber )

    if ( B % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( B % Type ), B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )

  end subroutine Finalize


end module Bundle_H__Form
