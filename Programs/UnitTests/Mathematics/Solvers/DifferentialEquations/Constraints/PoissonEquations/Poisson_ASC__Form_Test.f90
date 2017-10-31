module Poisson_ASC__Form_Test__Form

  use Basics
  use Manifolds
  use Poisson_ASC__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Poisson_ASC__Form_Test'

  type, public :: Poisson_ASC__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      Stream
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
!    type ( Poisson_ASC_Form ), allocatable :: &
!      Poisson
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Poisson_ASC__Form_Test_Form

contains


  subroutine Initialize ( PFT, Name )

    class ( Poisson_ASC__Form_Test_Form ) :: &
      PFT
    character ( * ), intent ( in ) :: &
      Name

  end subroutine Initialize


  subroutine Finalize ( PFT )

    type ( Poisson_ASC__Form_Test_Form ) :: &
      PFT

!    if ( allocated ( PFT % Laplacian ) ) &
!      deallocate ( PFT % Laplacian )
    if ( allocated ( PFT % Geometry ) ) &
      deallocate ( PFT % Geometry )
    if ( allocated ( PFT % Atlas ) ) &
      deallocate ( PFT % Atlas )
    if ( allocated ( PFT % Stream ) ) &
      deallocate ( PFT % Stream )

  end subroutine Finalize


end module Poisson_ASC__Form_Test__Form



program Poisson_ASC__Form_Test

  use Basics
  use Poisson_ASC__Form_Test__Form

  implicit none

  type ( Poisson_ASC__Form_Test_Form ), allocatable :: &
    PFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( PFT )
  call PFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( PFT )

  deallocate ( PROGRAM_HEADER )

end program Poisson_ASC__Form_Test
