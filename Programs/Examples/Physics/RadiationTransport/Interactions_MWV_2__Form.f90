module Interactions_MWV_2__Form

  !-- Interactions_MarshakWaveVaytet_2__Form, Vaytet et al. 2011
  
  use GenASiS
  use Interactions_MWV_1__Form

  implicit none
  private

  type, public, extends ( Interactions_MWV_1_Form ) :: Interactions_MWV_2_Form
    real ( KDR ) :: &
      SpecificOpacityFloor = 0.0_KDR, &
      EnergyMax = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, private, pass :: &
      Set_MWV_2_G
    procedure, private, pass :: &
      Set_MWV_2_S
    generic, public :: &
      Set => Set_MWV_2_G, Set_MWV_2_S
    final :: &
      Finalize
    procedure, private, pass :: &
      ComputeKernel_G
    procedure, private, pass :: &
      ComputeKernel_S
    procedure, private, pass :: &
      ComputeTimeScaleKernel_G
    procedure, private, pass :: &
      ComputeTimeScaleKernel_S
  end type Interactions_MWV_2_Form

    real ( KDR ), private, parameter :: &
      PlanckRatio  =  3.83223_KDR  !-- P_4 / P_3

contains


  subroutine InitializeAllocate_I &
               ( I, MomentsType, Units, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      MomentsType
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
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
      I % Type = 'an Interactions_MWV_2'

    call I % InitializeTemplate &
           ( MomentsType, Units, nValues, VariableOption, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_I


  subroutine Set_MWV_2_G ( I, Fluid, SpecificOpacity, EnergyMax )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      EnergyMax

    call I % Set ( Fluid, SpecificOpacity )

    I % EnergyMax = EnergyMax

  end subroutine Set_MWV_2_G


  subroutine Set_MWV_2_S &
               ( I, Fluid, Energy, d3_Energy, SpecificOpacity, &
                 SpecificOpacityFloor, EnergyMax, iBaseCell )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Energy, &
      d3_Energy
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity, &
      SpecificOpacityFloor, &
      EnergyMax
    integer ( KDI ), intent ( in ) :: &
      iBaseCell

    call I % Set ( Fluid, Energy, d3_Energy, SpecificOpacity, iBaseCell )

    I % SpecificOpacityFloor = SpecificOpacityFloor
    I % EnergyMax = EnergyMax

  end subroutine Set_MWV_2_S


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


  subroutine ComputeKernel_G ( I, TP, M, N, T, Xi_J, Chi_J, Chi_H, J_EQ )

    class ( Interactions_MWV_2_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      T
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H, &
      J_EQ

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      k_B, &
      a, &
      Kappa, &
      E_Max, &
      S, &
      S_EQ

    k_B    =  CONSTANT % BOLTZMANN
    a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa  =  I % SpecificOpacity
    E_Max  =  I % EnergyMax

    nValues  =  size ( Xi_J )

    !$OMP parallel do private ( iV, S, S_Eq ) 
    do iV = 1, nValues

      S     =  1.0_KDR  -  PlanckRatio * k_B * TP ( iV ) / E_Max
      S_EQ  =  1.0_KDR  -  PlanckRatio * k_B *  T ( iV ) / E_Max

      J_EQ  ( iV )  =  a  *  T ( iV ) ** 4

      Xi_J  ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  S_EQ  *  J_EQ ( iV )
      Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  S

      Chi_H ( iV )  =  Chi_J ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel_G


  subroutine ComputeKernel_S ( I, J_EQ, M, N, T, Xi_J, Chi_J, Chi_H )

    class ( Interactions_MWV_2_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_EQ
    real ( KDR ), intent ( in ) :: &
      M, &
      N, &
      T
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Kappa_0, &
      Kappa, &
      Kappa_Min, &
      E_Max

    associate ( E  =>  I % Energy )

    Kappa_0    =  I % SpecificOpacity
    Kappa_Min  =  I % SpecificOpacityFloor
    E_Max      =  I % EnergyMax
    nValues  =  size ( Xi_J )

    do iV = 1, nValues

      Kappa  =  max ( Kappa_0  *  ( 1.0_KDR  -  E ( iV ) / E_Max ), Kappa_Min )

      Xi_J  ( iV )  =  Kappa  *  M  *  N  *  J_EQ ( iV )
      Chi_J ( iV )  =  Kappa  *  M  *  N 
      Chi_H ( iV )  =  Chi_J ( iV )

    end do !-- iV

    end associate !-- E

  end subroutine ComputeKernel_S


  subroutine ComputeTimeScaleKernel_G ( I, TP, M, N, U, T, J, RT )

    class ( Interactions_MWV_2_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      U, &
      T, &
      J
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      RT

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      k_B, &
      a, &
      Kappa, &
      E_Max, &
      S, &
      S_EQ, &
      J_EQ, &
      Q, &
      SqrtTiny

    k_B    =  CONSTANT % BOLTZMANN
    a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa  =  I % SpecificOpacity
    E_Max  =  I % EnergyMax

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    nValues  =  size ( J )

    !$OMP parallel do private ( iV, S, S_EQ, J_EQ, Q ) 
    do iV = 1, nValues

      S     =  1.0_KDR  -  PlanckRatio * k_B * TP ( iV ) / E_Max
      S_EQ  =  1.0_KDR  -  PlanckRatio * k_B *  T ( iV ) / E_Max

      J_EQ  =  a  *  T ( iV ) ** 4

      Q     =  abs ( Kappa  *  M ( iV )  *  N ( iV )  &
                     *  ( S_EQ * J_EQ  -  S * J ( iV ) ) )

      RT ( iV )  =  U ( iV ) / max ( Q, SqrtTiny )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeTimeScaleKernel_G


  subroutine ComputeTimeScaleKernel_S ( I, J_EQ, J, dV, M, N, U, T, RT )

    class ( Interactions_MWV_2_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_EQ, &
      J, &
      dV
    real ( KDR ), intent ( in ) :: &
      M, &
      N, &
      U, &
      T
    real ( KDR ), intent ( out ) :: &
      RT

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Kappa_0, &
      Kappa, &
      Kappa_Min, &
      E_Max, &
      Q, &
      SqrtTiny

    associate ( E  =>  I % Energy )

    Kappa_0    =  I % SpecificOpacity
    Kappa_Min  =  I % SpecificOpacityFloor
    E_Max      =  I % EnergyMax

    SqrtTiny   =  sqrt ( tiny ( 0.0_KDR ) )

    nValues  =  size ( J )

    Q = 0.0_KDR
    do iV = 1, nValues

      Kappa  =  max ( Kappa_0  *  ( 1.0_KDR  -  E ( iV ) / E_Max ), Kappa_Min )

      Q  =  Q  +  abs ( Kappa * M * N * ( J_EQ ( iV ) - J ( iV ) ) ) &
                  *  dV ( iV )

    end do !-- iV
    RT  =  U / max ( Q, SqrtTiny )

    end associate !-- E

  end subroutine ComputeTimeScaleKernel_S


end module Interactions_MWV_2__Form
