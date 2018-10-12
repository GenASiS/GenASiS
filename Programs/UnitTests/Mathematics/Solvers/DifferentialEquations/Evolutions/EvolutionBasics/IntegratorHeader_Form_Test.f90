module IntegratorHeader_Form_Test__Form

  use Basics
  use IntegratorHeader_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'IntegratorHeader_Form_Test'

  type, public :: IntegratorHeader_Form_Test_Form
    type ( IntegratorHeaderForm ), allocatable :: &
      IntegratorHeader
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type IntegratorHeader_Form_Test_Form

contains

  
  subroutine Initialize ( IHFT )

    class ( IntegratorHeader_Form_Test_Form ), intent ( inout ), target :: &
      IHFT

    allocate ( IHFT % IntegratorHeader )
    associate ( I => IHFT % IntegratorHeader )
    call I % Initialize ( ProgramName )
    end associate !-- I

  end subroutine Initialize


  subroutine Finalize ( IHFT )

    type ( IntegratorHeader_Form_Test_Form ), intent ( inout ) :: &
      IHFT

    deallocate ( IHFT % IntegratorHeader )

  end subroutine Finalize


end module IntegratorHeader_Form_Test__Form



program IntegratorHeader_Form_Test

  use Basics
  use IntegratorHeader_Form_Test__Form
  
  implicit none

  type ( IntegratorHeader_Form_Test_Form ), allocatable :: &
    IHFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
  allocate ( IHFT )
  call IHFT % Initialize ( )
  deallocate ( IHFT )

  deallocate ( PROGRAM_HEADER )

end program IntegratorHeader_Form_Test
