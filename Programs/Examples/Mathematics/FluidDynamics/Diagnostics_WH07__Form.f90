module Diagnostics_WH07__Form

  use Basics

  implicit none
  private

  type, public, extends ( VariableGroupForm ) :: Diagnostics_WH07_Form
    integer ( KDI ) :: &
      N_FIELDS = 5, &
      PRESSURE_FORCE             = 1, &
      GRAVITATIONAL_FORCE_M      = 2, &
      GRAVITATIONAL_FORCE_PHI    = 3, &
      GRAVITATIONAL_FORCE_PHI_VJ = 4, &
      NET_FORCE                  = 5
  contains
    procedure, public, pass :: &
      Initialize_D
    generic, public :: &
      Initialize => Initialize_D
    procedure, public, pass :: &
      SetOutput
  end type Diagnostics_WH07_Form

contains


  subroutine Initialize_D ( D, Name, nValues )

    class ( Diagnostics_WH07_Form ), intent ( inout ) :: &
      D
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ) :: &
      nValues

    call D % VariableGroupForm % Initialize &
           ( [ nValues, D % N_FIELDS ], &
             VariableOption = [ 'PressureForce            ', &
                                'GravitationalForce_M     ', &
                                'GravitationalForce_Phi   ', &
                                'GravitationalForce_Phi_VJ', &
                                'NetForce                 ' ], &
             NameOption = Name, ClearOption = .true. )

  end subroutine Initialize_D


  subroutine SetOutput ( I, Output )

    class ( Diagnostics_WH07_Form ), intent ( inout ) :: &
      I
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize ( I )

  end subroutine SetOutput


end module Diagnostics_WH07__Form
