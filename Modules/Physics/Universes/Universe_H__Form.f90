module Universe_H__Form

  !-- Universe_Header__Form

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
    character ( LDL ), dimension ( : ), allocatable :: &
      dT_Label
    class ( Integrator_H_Form ), allocatable :: &
      Integrator
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, public, pass :: &
      Show => Show_U
    procedure, public, pass :: &
      Evolve
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass :: &
      ShowDiagnostics
  end type Universe_H_Form

    private :: &
      ShowSystem


contains

 
  subroutine Initialize_H ( U, Name )

    class ( Universe_H_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name

    U % IGNORABILITY = CONSOLE % INFO_1

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe' 

    U % Name  =  Name

    U % DeviceMemory  =  OffloadEnabled ( ) .and. NumberOfDevices ( ) >= 1
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

    call U % ShowParameters ( )
    call U % Integrator % Show ( )
    call U % ShowDiagnostics ( )

  end subroutine Show_U


  subroutine Evolve ( U )

    class ( Universe_H_Form ), intent ( inout ) :: &
      U

    associate ( I  =>  U % Integrator )
    if ( .not. associated ( I % ShowSystem ) ) &
      I % ShowSystem  =>  ShowSystem
    end associate !-- I

    call U % Integrator % Evolve ( )

  end subroutine Evolve


  subroutine Finalize ( U )

    type ( Universe_H_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Integrator ) ) &
      deallocate ( U % Integrator )
    if ( allocated ( U % dT_Label ) ) &
      deallocate ( U % dT_Label )

    if ( U % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( Universe_H_Form ), intent ( in ) :: &
      U

    call Show ( U % DeviceMemory, &
                'DeviceMemory', U % IGNORABILITY ) 
    call Show ( U % PinnedMemory, &
                'PinnedMemory', U % IGNORABILITY ) 
    call Show ( U % DevicesCommunicate, &
                'DevicesCommunicate', U % IGNORABILITY ) 

  end subroutine ShowParameters


  subroutine ShowDiagnostics ( U )

    class ( Universe_H_Form ), intent ( in ) :: &
      U

  end subroutine ShowDiagnostics


  subroutine ShowSystem ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    select type ( U  =>  I % System )
      class is ( Universe_H_Form )

    call U % Show ( )

    end select !-- U

  end subroutine ShowSystem


end module Universe_H__Form
