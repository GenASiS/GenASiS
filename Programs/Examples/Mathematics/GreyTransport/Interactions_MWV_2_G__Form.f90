module Interactions_MWV_2_G__Form

  !-- Interactions_MarshakWaveVaytet_2_Grey__Form, Vaytet et al. 2011
  
  use Basics
  use Fluid_P__Template
  use PhotonMoments_Form
  use Interactions_MWV_1_G__Form

  implicit none
  private

  type, public, extends ( Interactions_MWV_1_G_Form ) :: &
    Interactions_MWV_2_G_Form
      real ( KDR ) :: &
        EnergyMax = 0.0_KDR
  contains
    procedure, private, pass :: &
      Set_MWV_2_G
    generic, public :: &
      Set => Set_MWV_2_G
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeKernel
  end type Interactions_MWV_2_G_Form

    real ( KDR ), private, parameter :: &
      PlanckRatio  =  3.83223_KDR  !-- P_4 / P_3

contains


  subroutine InitializeAllocate_MWV_2_G &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_MWV_2_G_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_MWV_2_G'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, nValues, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_MWV_2_G


  subroutine Set_MWV_2_G ( I, Radiation, Fluid, SpecificOpacity, EnergyMax )

    class ( Interactions_MWV_2_G_Form ), intent ( inout ) :: &
      I
    class ( PhotonMomentsForm ), intent ( in ), target :: &
      Radiation
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      EnergyMax

    call I % Set ( Radiation, Fluid, SpecificOpacity )

    I % EnergyMax = EnergyMax

  end subroutine Set_MWV_2_G


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_2_G_Form ), intent ( inout ) :: &
      I

    !-- Trigger parent finalization
    
  end subroutine Finalize


  subroutine ComputeKernel ( TP, M, N, T, I, EDV, EOV, TOV )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      T
    class ( Interactions_MWV_2_G_Form ), intent ( in ) :: &
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
      S, &
      S_Eq

    k_B    =  CONSTANT % BOLTZMANN
    a      =  CONSTANT % RADIATION
    Kappa  =  I % SpecificOpacity
    E_Max  =  I % EnergyMax

    nValues  =  size ( EDV )

    !$OMP parallel do private ( iV, S, S_Eq ) 
    do iV = 1, nValues

      S     =  1.0_KDR  -  PlanckRatio * k_B * TP ( iV ) / E_Max
      S_Eq  =  1.0_KDR  -  PlanckRatio * k_B *  T ( iV ) / E_Max

      EDV ( iV )  =  a  *  T ( iV ) ** 4  *  S_Eq / S
      EOV ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  S
      TOV ( iV )  =  EOV ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_MWV_2_G__Form
