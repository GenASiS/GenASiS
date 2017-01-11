module AtlasHeader_Form_Test__Form

  use Basics
  use AtlasHeader_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'AtlasHeader_Form_Test'

  type, public :: AtlasHeader_Form_Test_Form
    type ( AtlasHeaderForm ), allocatable :: &
      AtlasHeaderBase, &
      AtlasHeaderFiber
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type AtlasHeader_Form_Test_Form

contains

  
  subroutine Initialize ( AHFT )

    class ( AtlasHeader_Form_Test_Form ), intent ( inout ), target :: &
      AHFT

    allocate ( AHFT % AtlasHeaderBase )
    associate ( A => AHFT % AtlasHeaderBase )
    call A % Initialize &
           ( 'Base', CommunicatorOption = PROGRAM_HEADER % Communicator, &
             iDimensionalityOption = 1 )
    call A % Connectivity % Show ( )
    call A % Show ( )
    end associate !-- AH

    allocate ( AHFT % AtlasHeaderFiber )
    associate ( A => AHFT % AtlasHeaderFiber )
    call A % Initialize &
           ( 'Fiber', iDimensionalityOption = 2 )
    call A % Connectivity % Show ( )
    call A % Show ( )
    end associate !-- AH

  end subroutine Initialize


  subroutine Finalize ( AHFT )

    type ( AtlasHeader_Form_Test_Form ), intent ( inout ) :: &
      AHFT

    deallocate ( AHFT % AtlasHeaderFiber )
    deallocate ( AHFT % AtlasHeaderBase )

  end subroutine Finalize


end module AtlasHeader_Form_Test__Form



program AtlasHeader_Form_Test

  use Basics
  use AtlasHeader_Form_Test__Form
  
  implicit none

  type ( AtlasHeader_Form_Test_Form ), allocatable :: &
    AHFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '3D_1D' )
    
  allocate ( AHFT )
  call AHFT % Initialize ( )
  deallocate ( AHFT )

  deallocate ( PROGRAM_HEADER )

end program AtlasHeader_Form_Test
