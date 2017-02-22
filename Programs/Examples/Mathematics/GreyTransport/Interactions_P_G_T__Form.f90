module Interactions_P_G_T__Form

  !-- Interactions_Photons_Grey_Temperature__Form
  
  use Basics
  use Fluid_P__Template
  use PhotonMoments_Form
  use Interactions_P_G_L__Form

  implicit none
  private

  type, public, extends ( Interactions_P_G_L_Form ) :: Interactions_P_G_T_Form
    real ( KDR ) :: &
      TemperatureScale = 1.0_KDR
  contains
    procedure, private, pass :: &
      Set_P_G_T
    generic, public :: &
      Set => Set_P_G_T
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeKernel
  end type Interactions_P_G_T_Form

    real ( KDR ), private, parameter :: &
      PlanckRatio  =  3.83223_KDR  !-- P_4 / P_3

contains


  subroutine InitializeAllocate_P_G_T &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_P_G_T_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_P_G_T'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, nValues, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_P_G_T


  subroutine Set_P_G_T &
               ( I, Radiation, Fluid, SpecificOpacity, EnergyMax, &
                 TemperatureScale )

    class ( Interactions_P_G_T_Form ), intent ( inout ) :: &
      I
    class ( PhotonMomentsForm ), intent ( in ), target :: &
      Radiation
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      EnergyMax, &
      TemperatureScale

    call I % Set ( Radiation, Fluid, SpecificOpacity, EnergyMax )

    I % TemperatureScale = TemperatureScale

  end subroutine Set_P_G_T


  impure elemental subroutine Finalize ( I )

    type ( Interactions_P_G_T_Form ), intent ( inout ) :: &
      I

    !-- Trigger parent finalization
    
  end subroutine Finalize


  subroutine ComputeKernel ( TP, M, N, T, I, EDV, EOV, TOV )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      T
    class ( Interactions_P_G_T_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EDV, &
      EOV, &
      TOV

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      k_B, &
      a, &
      Kappa, &
      E_Max, &
      T_0, &
      S, &
      S_Eq

    k_B    =  CONSTANT % BOLTZMANN
    a      =  CONSTANT % RADIATION
    Kappa  =  I % SpecificOpacity
    E_Max  =  I % EnergyMax
    T_0    =  I % TemperatureScale

    nValues  =  size ( EDV )

    !$OMP parallel do private ( iV, S, S_Eq ) 
    do iV = 1, nValues

      S     =  1.0_KDR  -  PlanckRatio * k_B * TP ( iV ) / E_Max
      S_Eq  =  1.0_KDR  -  PlanckRatio * k_B *  T ( iV ) / E_Max

      EDV ( iV )  =  a  *  T ( iV ) ** 4  *  S_Eq / S
      EOV ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  &
                     *  ( T ( iV ) / T_0 ) ** 1.5_KDR  *  S
      TOV ( iV )  =  EOV ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_P_G_T__Form
