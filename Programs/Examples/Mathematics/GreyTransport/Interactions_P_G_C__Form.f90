module Interactions_P_G_C__Form

  !-- Interactions_Photons_Grey_Constant__Form
  
  use Basics
  use Fluid_P__Template
  use PhotonMoments_Form
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_P_G_C_Form
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR
    class ( Fluid_P_Template ), pointer :: &
      Fluid => null ( )
    class ( PhotonMomentsForm ), pointer :: &
      Radiation => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_P_G_C
    generic, public :: &
      Initialize => InitializeAllocate_P_G_C
    procedure, private, pass :: &
      Set_P_G_C
    generic, public :: &
      Set => Set_P_G_C
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeKernel
  end type Interactions_P_G_C_Form

contains


  subroutine InitializeAllocate_P_G_C &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_P_G_C_Form ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_P_G_C'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, nValues, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_P_G_C


  subroutine Set_P_G_C ( I, Radiation, Fluid, SpecificOpacity )

    class ( Interactions_P_G_C_Form ), intent ( inout ) :: &
      I
    class ( PhotonMomentsForm ), intent ( in ), target :: &
      Radiation
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

    I % SpecificOpacity  =   SpecificOpacity
    I % Fluid            =>  Fluid
    I % Radiation        =>  Radiation

  end subroutine Set_P_G_C


  subroutine Compute ( I )

    class ( Interactions_P_G_C_Form ), intent ( inout ) :: &
      I

    associate &
      ( R => I % Radiation, &
        F => I % Fluid )
    call I % ComputeKernel &
           ( R % Value ( :, R % TEMPERATURE_PARAMETER ), &
             F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % COMOVING_DENSITY ), &
             F % Value ( :, F % TEMPERATURE ), &
             I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
             I % Value ( :, I % EFFECTIVE_OPACITY ), &
             I % Value ( :, I % TRANSPORT_OPACITY ) )
    end associate !-- R, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_P_G_C_Form ), intent ( inout ) :: &
      I

    nullify ( I % Fluid )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( TP, M, N, T, I, EDV, EOV, TOV )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      T
    class ( Interactions_P_G_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EDV, &
      EOV, &
      TOV

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a, &
      Kappa

    a      =  CONSTANT % RADIATION
    Kappa  =  I % SpecificOpacity

    nValues  =  size ( EDV )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      EDV ( iV )  =  a  *  T ( iV ) ** 4
      EOV ( iV )  =  Kappa  *  M ( iV )  *  N ( iV ) 
      TOV ( iV )  =  EOV ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_P_G_C__Form
