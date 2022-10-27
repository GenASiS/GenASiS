module Universe_R_B__Form

  !-- Universe_Radiation_Box__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Universe_F_B__Form

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: Universe_R_B_Form
    character ( LDL ) :: &
      MomentsType = ''
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
  contains
    procedure, private, pass :: &
      Initialize_R_B
    generic, public :: &
      Initialize => Initialize_R_B
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type Universe_R_B_Form


contains


  subroutine Initialize_R_B &
               ( U, RadiationName, RadiationType, MomentsType, Name )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      MomentsType, &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_B'

    call U % Universe_H_Form % Initialize ( Name )

    U % MomentsType  =  MomentsType

    allocate ( U % Units_R ( 1 ) )
    call U % Units_R ( 1 ) % Initialize ( )

  end subroutine Initialize_R_B


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( Universe_R_B_Form ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % MomentsType, 'MomentsType', U % IGNORABILITY )

  end subroutine ShowParameters


end module Universe_R_B__Form
