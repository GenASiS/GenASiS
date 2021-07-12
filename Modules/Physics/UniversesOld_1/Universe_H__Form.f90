module Universe_H__Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public :: Universe_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    logical ( KDL ) :: &
      DeviceMemory, &
      PinnedMemory, &
      DevicesCommunicate
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    class ( Integrator_H_Form ), allocatable :: &
      Integrator
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, private, pass :: &
      Show_U
    generic, public :: &
      Show => Show_U
    procedure, public, pass :: &
      Evolve
    final :: &
      Finalize
  end type Universe_H_Form

contains

 
  subroutine Initialize_H ( U, NameOption )

    class ( Universe_H_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    U % IGNORABILITY = CONSOLE % INFO_1

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe' 

    U % Name  =  'Universe'
    if ( present ( NameOption ) ) &
      U % Name  =  NameOption

    U % DeviceMemory  =  OffloadEnabled ( ) .and. GetNumberOfDevices ( ) >= 1
    call PROGRAM_HEADER % GetParameter ( U % DeviceMemory, 'DeviceMemory' )

    U % PinnedMemory        =  U % DeviceMemory
    U % DevicesCommunicate  =  U % DeviceMemory
    call PROGRAM_HEADER % GetParameter &
           ( U % PinnedMemory, 'PinnedMemory' )
    call PROGRAM_HEADER % GetParameter &
           ( U % DevicesCommunicate, 'DevicesCommunicate' )

    call Show ( 'Initializing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

  end subroutine Initialize_H


  subroutine Show_U ( U )

    class ( Universe_H_Form ), intent ( in ) :: &
      U

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( U % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

    call Show ( U % DeviceMemory, &
                'DeviceMemory', U % IGNORABILITY ) 
    call Show ( U % PinnedMemory, &
                'PinnedMemory', U % IGNORABILITY ) 
    call Show ( U % DevicesCommunicate, &
                'DevicesCommunicate', U % IGNORABILITY ) 

    call U % Integrator % Show ( )

  end subroutine Show_U


  subroutine Evolve ( U )

    class ( Universe_H_Form ), intent ( inout ) :: &
      U

    call U % Integrator % Evolve ( )

  end subroutine Evolve


  subroutine Finalize ( U )

    type ( Universe_H_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Integrator ) ) &
      deallocate ( U % Integrator )

    if ( U % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

  end subroutine Finalize


end module Universe_H__Form
