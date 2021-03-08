module Atlas_SC_SC__Form_Test__Form

  use Basics
  use Charts
  use Atlas_SC_SC__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Atlas_SC_SC__Form_Test'

  type, public :: Atlas_SC_SC_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GIS
    type ( Atlas_SC_SC_Form ), allocatable :: &
      Atlas
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Atlas_SC_SC_Form_Test_Form

contains

  
  subroutine Initialize ( AFT, Name )

    class ( Atlas_SC_SC_Form_Test_Form ), intent ( inout ), target :: &
      AFT
    character ( * ), intent ( in ) :: &
      Name

    allocate &
      ( AFT % Atlas, &
        AFT % GIS )
    associate &
      ( A   => AFT % Atlas, &
        GIS => AFT % GIS )

    call A % Initialize &
           ( 'Atlas', CommunicatorOption = PROGRAM_HEADER % Communicator )
    call A % CreateChart_SC ( )  
    call A % SetGeometry ( )

    call GIS % Initialize &
           ( PROGRAM_HEADER % Name, CommunicatorOption = A % Communicator )
    call A % OpenStream ( GIS, '1', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    end associate !-- A, etc.

  end subroutine Initialize


  subroutine Finalize ( AFT )

    type ( Atlas_SC_SC_Form_Test_Form ), intent ( inout ) :: &
      AFT

    deallocate ( AFT % Atlas )
    deallocate ( AFT % GIS )

  end subroutine Finalize


end module Atlas_SC_SC__Form_Test__Form



program Atlas_SC_SC__Form_Test

  use Basics
  use Atlas_SC_SC__Form_Test__Form
  
  implicit none

  type ( Atlas_SC_SC_Form_Test_Form ), allocatable :: &
    AFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '2D' )

  allocate ( AFT )
  call AFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( AFT )

  deallocate ( PROGRAM_HEADER )

end program Atlas_SC_SC__Form_Test
