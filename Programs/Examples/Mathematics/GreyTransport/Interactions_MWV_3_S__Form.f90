module Interactions_MWV_3_S__Form

  !-- Interactions_MarshakWaveVaytet_3_Spectral__Form, Vaytet et al. 2011
  
  use Basics
  use Mathematics
  use Fluid_P__Template
  use SetPlanckSpectrum_Command
  use PhotonMoments_S__Form
  use Interactions_MWV_2_S__Form

  implicit none
  private

  type, public, extends ( Interactions_MWV_2_S_Form ) :: &
    Interactions_MWV_3_S_Form
      real ( KDR ) :: &
        TemperatureScale = 1.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_MWV_1_S
    procedure, private, pass :: &
      Set_MWV_3_S
    generic, public :: &
      Set => Set_MWV_3_S
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeKernel
  end type Interactions_MWV_3_S_Form

contains


  subroutine InitializeAllocate_MWV_1_S &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_MWV_3_S_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_MWV_3_S'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_MWV_1_S


  subroutine Set_MWV_3_S &
               ( I, Radiation, Fluid, Energy, SpecificOpacity, EnergyMax, &
                 TemperatureScale, iBaseCell )

    class ( Interactions_MWV_3_S_Form ), intent ( inout ) :: &
      I
    class ( PhotonMoments_S_Form ), intent ( in ), target :: &
      Radiation
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Energy
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      EnergyMax, &
      TemperatureScale
    integer ( KDI ), intent ( in ) :: &
      iBaseCell

    call I % Set ( Radiation, Fluid, Energy, SpecificOpacity, EnergyMax, &
                   iBaseCell )

    I % TemperatureScale = TemperatureScale

  end subroutine Set_MWV_3_S


  subroutine Compute ( I, Current )

    class ( Interactions_MWV_3_S_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    associate &
      ( R   => I % Radiation, &
        F   => I % Fluid, &
        iBC => I % iBaseCell )
    call I % ComputeKernel &
           ( F % Value ( iBC, F % BARYON_MASS ), &
             F % Value ( iBC, F % COMOVING_DENSITY ), &
             F % Value ( iBC, F % TEMPERATURE ), &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ), &
             R % Value ( :, R % COMOVING_ENERGY_EQ ) )
    end associate !-- R, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_3_S_Form ), intent ( inout ) :: &
      I

    nullify ( I % Radiation )
    nullify ( I % Fluid )
    nullify ( I % Energy )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( M, N, T, I, Xi_J, Chi_J, Chi_H, J_Eq )

    real ( KDR ), intent ( in ) :: &
      M, &
      N, &
      T
    class ( Interactions_MWV_3_S_Form ), intent ( in ) :: &
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
      Kappa_0, &
      Kappa, &
      Kappa_Min, &
      E_Max, &
      T_0

    associate ( E  =>  I % Energy )

    call SetPlanckSpectrum ( E, T, J_Eq )

    Kappa_0    =  I % SpecificOpacity
    Kappa_Min  =  10.0_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    E_Max      =  I % EnergyMax
    T_0        =  I % TemperatureScale
    nValues    =  size ( Xi_J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      Kappa  =  max ( Kappa_0  *  ( 1.0_KDR  -  E ( iV ) / E_Max ), Kappa_Min )

      Xi_J  ( iV )  =  Kappa  *  M  *  N  *  ( T / T_0 ) ** 1.5_KDR  &
                       *  J_Eq ( iV )
      Chi_J ( iV )  =  Kappa  *  M  *  N  *  ( T / T_0 ) ** 1.5_KDR
      Chi_H ( iV )  =  Chi_J ( iV )

    end do !-- iV
    !$OMP end parallel do

    end associate !-- E

  end subroutine ComputeKernel


end module Interactions_MWV_3_S__Form
