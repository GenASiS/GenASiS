module Diagnostics_WH07__Form

  use Basics

  implicit none
  private

  type, public, extends ( VariableGroupForm ) :: Diagnostics_WH07_Form
    integer ( KDI ) :: &
      N_FIELDS = 8, &
      PRESSURE_FORCE             = 1, &
      GRAVITATIONAL_POTENTIAL_C  = 2, &
      GRAVITATIONAL_FORCE_M      = 3, &
      GRAVITATIONAL_FORCE_PHI    = 4, &
      GRAVITATIONAL_FORCE_PHI_VJ = 5, &
      GRAVITATIONAL_FORCE_PHI_C  = 6, &
      NET_FORCE                  = 7, &
      NET_POWER                  = 8
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
                                'GravitationalPotential_C ', &
                                'GravitationalForce_M     ', &
                                'GravitationalForce_Phi   ', &
                                'GravitationalForce_Phi_VJ', &
                                'GravitationalForce_Phi_C ', &
                                'NetForce                 ', &
                                'NetPower                 ' ], &
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
