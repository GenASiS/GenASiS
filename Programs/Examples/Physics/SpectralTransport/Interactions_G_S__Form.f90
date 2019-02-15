module Interactions_G_S__Form

  !-- Interactions_Generic_Spectral_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_G_S_Form
    integer ( KDI ) :: &
      iBaseCell = 0
    real ( KDR ) :: &
      OpacityAbsorption = 0.0_KDR
    real ( KDR ), dimension ( : ), pointer :: &
      Energy => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      Set => Set_G_S
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private, pass :: &
      ComputeKernel
  end type Interactions_G_S_Form

contains


  subroutine InitializeAllocate_I &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_G_S_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_G_S'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_I


  subroutine Set_G_S ( I, Fluid, Energy, OpacityAbsorption, iBaseCell )

    class ( Interactions_G_S_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Energy
    real ( KDR ), intent ( in ) :: &
      OpacityAbsorption
    integer ( KDI ), intent ( in ) :: &
      iBaseCell

    I % iBaseCell          =   iBaseCell
    I % OpacityAbsorption  =   OpacityAbsorption
    I % Energy             =>  Energy
    I % Fluid              =>  Fluid

  end subroutine Set_G_S


  subroutine Compute ( I, R )

    class ( Interactions_G_S_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( inout ) :: &
      R

    associate &
      (   F => I % Fluid, &
        iBC => I % iBaseCell )
    call SetPlanckSpectrum &
           ( I % Energy, &
             F % Value ( iBC, F % TEMPERATURE ), &
             I % Value ( :, I % EQUILIBRIUM_J ) )
    end associate !-- F, etc.

    call I % ComputeKernel &
           ( I % Value ( :, I % EQUILIBRIUM_J ), &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ) )

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_G_S_Form ), intent ( inout ) :: &
      I

    nullify ( I % Energy )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( I, Xi_J, Chi_J, Chi_H, J_Eq )

    class ( Interactions_G_S_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_Eq
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Kappa_A

    Kappa_A  =  I % OpacityAbsorption
    nValues  =  size ( Xi_J )

    do iV = 1, nValues
      Xi_J  ( iV )  =  Kappa_A  *  J_Eq ( iV )
      Chi_J ( iV )  =  Kappa_A
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV

  end subroutine ComputeKernel


end module Interactions_G_S__Form
