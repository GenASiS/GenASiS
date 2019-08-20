module Interactions_C__Form

  !-- Interactions_Constant_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_C_Form
    real ( KDR ) :: &
      OpacityAbsorption = 0.0_KDR
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
  end type Interactions_C_Form

contains


  subroutine InitializeAllocate_I &
               ( I, MomentsType, Units, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_C_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_C'

    call I % InitializeTemplate &
           ( MomentsType, Units, nValues, VariableOption, NameOption, &
             ClearOption, UnitOption )

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
    select type ( SF => F % Sources )
    class is ( Sources_F_Form )

    select type ( R )
    class is ( RadiationMomentsForm )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      call I % ComputeTimeScaleKernel_G &
             (  F % Value ( :,  F % INTERNAL_ENERGY ), &
                F % Value ( :,  F % TEMPERATURE ), &
                R % Value ( :,  R % COMOVING_ENERGY ), &
               SF % Value ( :, SF % RADIATION_TIME ) )
    case ( 'SPECTRAL' )
      call SetPlanckSpectrum &
             ( I % Energy, &
               F % Value ( iBC, F % TEMPERATURE ), &
               I % Value ( :, I % EQUILIBRIUM_J ) )
      call I % ComputeTimeScaleKernel_S &
             (  I % Value ( :, I % EQUILIBRIUM_J ), &
                R % Value ( :, R % COMOVING_ENERGY ), &
                I % d3_Energy, &
                F % Value ( iBC,  F % INTERNAL_ENERGY ), &
               SF % Value ( iBC, SF % RADIATION_TIME ) )
    end select !-- MomentsType

    end select !-- R
    end select !-- SF
    end associate !-- F, etc.

  end subroutine ComputeTimeScale


  impure elemental subroutine Finalize ( I )

    type ( Interactions_C_Form ), intent ( inout ) :: &
      I

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel_G ( I, T, Xi_J, Chi_J, Chi_H, J_EQ )

    class ( Interactions_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
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
      a, &
      Kappa_A

    a        =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa_A  =  I % OpacityAbsorption

    nValues  =  size ( Xi_J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      J_EQ  ( iV )  =  a  *  T ( iV ) ** 4
      Xi_J  ( iV )  =  Kappa_A  *  J_EQ ( iV )
      Chi_J ( iV )  =  Kappa_A
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel_G


  subroutine ComputeKernel_S ( I, J_EQ, Xi_J, Chi_J, Chi_H )

    class ( Interactions_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_EQ
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
      Xi_J  ( iV )  =  Kappa_A  *  J_EQ ( iV )
      Chi_J ( iV )  =  Kappa_A
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV

  end subroutine ComputeKernel_S


  subroutine ComputeTimeScaleKernel_G ( I, U, T, J, RT )

    class ( Interactions_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      U, &
      T, &
      J
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      RT

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a, &
      Kappa_A, &
      J_EQ, &
      Q, &
      SqrtTiny

    a        =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa_A  =  I % OpacityAbsorption

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    nValues  =  size ( J )

    !$OMP parallel do private ( iV, J_EQ, Q ) 
    do iV = 1, nValues
      J_EQ       =  a  *  T ( iV ) ** 4
      Q          =  abs ( Kappa_A * ( J_EQ - J ( iV ) ) )
      RT ( iV )  =  U ( iV ) / max ( Q, SqrtTiny )
    end do !-- iV
    !$OMP end parallel do


  end subroutine ComputeTimeScaleKernel_G


  subroutine ComputeTimeScaleKernel_S ( I, J_EQ, J, dV, U, RT )

    class ( Interactions_C_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_EQ, &
      J, &
      dV
    real ( KDR ), intent ( in ) :: &
      U
    real ( KDR ), intent ( out ) :: &
      RT

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Kappa_A, &
      Q, &
      SqrtTiny

    Kappa_A  =  I % OpacityAbsorption

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    nValues  =  size ( J )

    Q = 0.0_KDR
    do iV = 1, nValues
      Q  =  Q  +  abs ( Kappa_A * ( J_EQ ( iV ) - J ( iV ) ) )  *  dV ( iV )
    end do !-- iV
    RT  =  U / max ( Q, SqrtTiny )

  end subroutine ComputeTimeScaleKernel_S


end module Interactions_C__Form
