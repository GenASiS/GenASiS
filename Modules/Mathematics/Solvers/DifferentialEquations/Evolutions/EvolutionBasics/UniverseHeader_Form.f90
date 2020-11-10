!-- UniverseHeader handles metadata of a Universe.

module UniverseHeader_Form

  use Basics

  implicit none
  private

  type, public :: UniverseHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0
    logical ( KDL ) :: &
      UseDevice = .false., &
      UseCustomBoundaryInner = .false.
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( Character_1D_Form ), dimension ( : ), allocatable :: &
      BoundaryConditionsFace
    procedure ( SI ), pointer :: &
      SetInitial => null ( )
    procedure ( RI ), pointer :: &
      ResetInitial => null ( )
  contains
    procedure, public, pass :: &
      InitializeHeader
    final :: &
      Finalize
  end type UniverseHeaderForm

  
  interface

    subroutine SI ( U )
      import UniverseHeaderForm
      class ( UniverseHeaderForm ), intent ( inout ) :: &
        U
    end subroutine SI

    subroutine RI ( U, RestartFrom, RestartTime )
      use Basics
      import UniverseHeaderForm
      class ( UniverseHeaderForm ), intent ( inout ) :: &
        U
      integer ( KDI ), intent ( in ) :: &
        RestartFrom
      real ( KDR ), intent ( out ) :: &
        RestartTime
    end subroutine RI

  end interface


contains


  subroutine InitializeHeader ( U, Name )

    class ( UniverseHeaderForm ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      Name

    U % IGNORABILITY = CONSOLE % INFO_1

    if ( U % Type == '' ) &
      U % Type = 'a Universe' 

    U % Name = Name

    U % UseDevice = ( OffloadEnabled ( ) .and. GetNumberOfDevices ( ) >= 1 )
    call PROGRAM_HEADER % GetParameter ( U % UseDevice, 'UseDevice' )

    call Show ( 'Initializing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )
    call Show ( U % UseDevice, 'UseDevice', U % IGNORABILITY )

  end subroutine InitializeHeader


  subroutine Finalize ( U )

    type ( UniverseHeaderForm ), intent ( inout ) :: &
      U

    if ( allocated ( U % BoundaryConditionsFace ) ) &
      deallocate ( U % BoundaryConditionsFace )

    if ( U % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( U % Type ), U % IGNORABILITY )
    call Show ( U % Name, 'Name', U % IGNORABILITY )

  end subroutine Finalize


end module UniverseHeader_Form
