module Interactions_MWV_1_S__Form

  !-- Interactions_MarshakWaveVaytet_1_Spectral__Form, Vaytet et al. 2011
  
  use Basics
  use Mathematics
  use Fluid_P__Template
  use SetPlanckSpectrum_Command
  use PhotonMoments_S__Form
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_MWV_1_S_Form
    integer ( KDI ) :: &
      iBaseCell = 0
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR
    real ( KDR ), dimension ( : ), pointer :: &
      Energy => null ( )
    class ( Fluid_P_Template ), pointer :: &
      Fluid => null ( )
    class ( PhotonMoments_S_Form ), pointer :: &
      Radiation => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_MWV_1_S
    generic, public :: &
      Initialize => InitializeAllocate_MWV_1_S
    procedure, private, pass :: &
      Set_MWV_1_S
    generic, public :: &
      Set => Set_MWV_1_S
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeKernel
  end type Interactions_MWV_1_S_Form

contains


  subroutine InitializeAllocate_MWV_1_S &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_MWV_1_S_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_MWV_1_S'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_MWV_1_S


  subroutine Set_MWV_1_S &
               ( I, Radiation, Fluid, Energy, SpecificOpacity, iBaseCell )

    class ( Interactions_MWV_1_S_Form ), intent ( inout ) :: &
      I
    class ( PhotonMoments_S_Form ), intent ( in ), target :: &
      Radiation
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Energy
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity
    integer ( KDI ), intent ( in ) :: &
      iBaseCell

    I % iBaseCell        =   iBaseCell
    I % SpecificOpacity  =   SpecificOpacity
    I % Energy           =>  Energy
    I % Fluid            =>  Fluid
    I % Radiation        =>  Radiation

  end subroutine Set_MWV_1_S


  subroutine Compute ( I, Current )

    class ( Interactions_MWV_1_S_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    associate &
      ( R   => I % Radiation, &
        F   => I % Fluid, &
        iBC => I % iBaseCell )

    call SetPlanckSpectrum &
           ( I % Energy, &
             F % Value ( iBC, F % TEMPERATURE ), &
             R % Value ( :, R % COMOVING_ENERGY_EQ ) )

    call I % ComputeKernel &
           ( R % Value ( :, R % COMOVING_ENERGY_EQ ), &
             F % Value ( iBC, F % BARYON_MASS ), &
             F % Value ( iBC, F % COMOVING_DENSITY ), &
             F % Value ( iBC, F % TEMPERATURE ), &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ) )

    end associate !-- R, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_1_S_Form ), intent ( inout ) :: &
      I

    nullify ( I % Radiation )
    nullify ( I % Fluid )
    nullify ( I % Energy )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( J_Eq, M, N, T, I, Xi_J, Chi_J, Chi_H )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_Eq
    real ( KDR ), intent ( in ) :: &
      M, &
      N, &
      T
    class ( Interactions_MWV_1_S_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Kappa

    Kappa    =  I % SpecificOpacity
    nValues  =  size ( Xi_J )

    do iV = 1, nValues
      Xi_J  ( iV )  =  Kappa  *  M  *  N  *  J_Eq ( iV )
      Chi_J ( iV )  =  Kappa  *  M  *  N 
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV

  end subroutine ComputeKernel


end module Interactions_MWV_1_S__Form
