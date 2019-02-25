module Interactions_C__Form

  !-- Interactions_Constant_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_C_Form
    integer ( KDI ) :: &
      iBaseCell = 0
    real ( KDR ) :: &
      OpacityAbsorption = 0.0_KDR
    real ( KDR ), dimension ( : ), pointer :: &
      Energy => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, private, pass :: &
      Set_C_G
    procedure, private, pass :: &
      Set_C_S
    generic, public :: &
      Set => Set_C_G, Set_C_S
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      ComputeTimeScale
    procedure, private, pass :: &
      ComputeKernel_G
    procedure, private, pass :: &
      ComputeKernel_S
    procedure, private, pass :: &
      ComputeTimeScaleKernel_G
  end type Interactions_C_Form

contains


  subroutine InitializeAllocate_I &
               ( I, MomentsType, LengthUnit, EnergyDensityUnit, &
                 TemperatureUnit, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    character ( * ), intent ( in ) :: &
      MomentsType
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
      I % Type = 'an Interactions_C'

    call I % InitializeTemplate &
           ( MomentsType, LengthUnit, EnergyDensityUnit, TemperatureUnit, &
             nValues, VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_I


  subroutine Set_C_G ( I, Fluid, OpacityAbsorption )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      OpacityAbsorption

    I % OpacityAbsorption  =   OpacityAbsorption
    I % Fluid              =>  Fluid

  end subroutine Set_C_G


  subroutine Set_C_S ( I, Fluid, Energy, OpacityAbsorption, iBaseCell )

    class ( Interactions_C_Form ), intent ( inout ) :: &
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

  end subroutine Set_C_S


  subroutine Compute ( I, R )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    associate &
      (   F => I % Fluid, &
        iBC => I % iBaseCell )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      call I % ComputeKernel_G &
             ( F % Value ( :, F % TEMPERATURE ), &
               I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ), &
               I % Value ( :, I % EQUILIBRIUM_J ) )
    case ( 'SPECTRAL' )
      call SetPlanckSpectrum &
             ( I % Energy, &
               F % Value ( iBC, F % TEMPERATURE ), &
               I % Value ( :, I % EQUILIBRIUM_J ) )
      call I % ComputeKernel_S &
             ( I % Value ( :, I % EQUILIBRIUM_J ), &
               I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ) )
    end select !-- MomentsType

    end associate !-- F, etc.

  end subroutine Compute


  subroutine ComputeTimeScale ( I, R )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    associate &
      (   F => I % Fluid, &
        iBC => I % iBaseCell )

    select type ( R )
    class is ( RadiationMomentsForm )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      call I % ComputeTimeScaleKernel_G &
             ( F % Value ( :, F % TEMPERATURE ), &
               R % Value ( :, R % COMOVING_ENERGY ), &
               I % Value ( :, I % TIME_SCALE ) )
    case ( 'SPECTRAL' )
    end select !-- MomentsType

    end select !-- R
    end associate !-- F, etc.

  end subroutine ComputeTimeScale


  subroutine ComputeKernel_G ( I, T, Xi_J, Chi_J, Chi_H, J_Eq )

    class ( Interactions_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      T
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
      Kappa_A

    a        =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa_A  =  I % OpacityAbsorption

    nValues  =  size ( Xi_J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
      Xi_J  ( iV )  =  Kappa_A  *  J_Eq ( iV )
      Chi_J ( iV )  =  Kappa_A
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel_G


  subroutine ComputeKernel_S ( I, J_Eq, Xi_J, Chi_J, Chi_H )

    class ( Interactions_C_Form ), intent ( in ) :: &
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

  end subroutine ComputeKernel_S


  subroutine ComputeTimeScaleKernel_G ( I, T, J, TS )

    class ( Interactions_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      T, &
      J
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      TS

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a, &
      Kappa_A, &
      J_Eq, &
      Q

    a        =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa_A  =  I % OpacityAbsorption

    nValues  =  size ( TS )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      J_Eq       =  a  *  T ( iV ) ** 4
      Q          =  abs ( Kappa_A * ( J_Eq - J ( iV ) ) )
      TS ( iV )  =  J_Eq / Q
    end do !-- iV
    !$OMP end parallel do


  end subroutine ComputeTimeScaleKernel_G


end module Interactions_C__Form
