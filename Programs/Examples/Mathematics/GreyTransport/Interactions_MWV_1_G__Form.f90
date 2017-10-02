module Interactions_MWV_1_G__Form

  !-- Interactions_MarshakWaveVaytet_1_Grey__Form, Vaytet et al. 2011
  
  use Basics
  use Mathematics
  use Fluid_P__Template
  use PhotonMoments_G__Form
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_MWV_1_G_Form
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR
    class ( Fluid_P_Template ), pointer :: &
      Fluid => null ( )
    class ( PhotonMoments_G_Form ), pointer :: &
      Radiation => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_MWV_1_G
    generic, public :: &
      Initialize => InitializeAllocate_MWV_1_G
    procedure, private, pass :: &
      Set_MWV_1_G
    generic, public :: &
      Set => Set_MWV_1_G
    procedure, public, pass :: &
      Compute
    procedure, private, pass ( I ) :: &
      ComputeEquilibrium_T
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeKernel
  end type Interactions_MWV_1_G_Form

contains


  subroutine InitializeAllocate_MWV_1_G &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_MWV_1_G_Form ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_MWV_1_G'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_MWV_1_G


  subroutine Set_MWV_1_G ( I, Radiation, Fluid, SpecificOpacity )

    class ( Interactions_MWV_1_G_Form ), intent ( inout ) :: &
      I
    class ( PhotonMoments_G_Form ), intent ( in ), target :: &
      Radiation
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

    I % SpecificOpacity  =   SpecificOpacity
    I % Fluid            =>  Fluid
    I % Radiation        =>  Radiation

  end subroutine Set_MWV_1_G


  subroutine Compute ( I, Current )

    class ( Interactions_MWV_1_G_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    associate &
      ( R => I % Radiation, &
        F => I % Fluid )
    call I % ComputeKernel &
           ( R % Value ( :, R % TEMPERATURE_PARAMETER ), &
             F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % COMOVING_DENSITY ), &
             F % Value ( :, F % TEMPERATURE ), &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ), &
             I % Value ( :, I % EQUILIBRIUM_J ) )
    end associate !-- R, etc.

  end subroutine Compute


  subroutine ComputeEquilibrium_T ( T_EQ, I, C )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T_EQ
    class ( Interactions_MWV_1_G_Form ), intent ( in ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      C

    if ( .not.associated ( I % Radiation ) &
         .or. .not.associated ( I % Fluid ) ) &
      return

    associate &
      ( R => I % Radiation, &
        F => I % Fluid )
    call Copy ( F % Value ( :, F % TEMPERATURE ), &
                R % Value ( :, R % TEMPERATURE_PARAMETER_EQ ) )
    end associate !-- R, etc.

  end subroutine ComputeEquilibrium_T


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_1_G_Form ), intent ( inout ) :: &
      I

    nullify ( I % Radiation )
    nullify ( I % Fluid )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( TP, M, N, T, I, Xi_J, Chi_J, Chi_H, J_Eq )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      T
    class ( Interactions_MWV_1_G_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H, &
      J_Eq

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a, &
      Kappa

    a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa  =  I % SpecificOpacity

    nValues  =  size ( Xi_J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      J_Eq  ( iV )  =  a  *  T ( iV ) ** 4

      Xi_J  ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  J_Eq ( iV )
      Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV ) 

      Chi_H ( iV )  =  Chi_J ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_MWV_1_G__Form
