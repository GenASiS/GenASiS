module Interactions_MWV_1__Form

  !-- Interactions_MarshakWaveVaytet_1__Form, Vaytet et al. 2011
  
  use GenASiS

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_MWV_1_Form
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, private, pass :: &
      Set_MWV_1_G
    procedure, private, pass :: &
      Set_MWV_1_S
    generic, public :: &
      Set => Set_MWV_1_G, Set_MWV_1_S
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
    procedure, private, pass :: &
      ComputeTimeScaleKernel_S
  end type Interactions_MWV_1_Form

contains


  subroutine InitializeAllocate_I &
               ( I, MomentsType, LengthUnit, EnergyDensityUnit, &
                 TemperatureUnit, nValues, VariableOption, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_MWV_1'

    call I % InitializeTemplate &
           ( MomentsType, LengthUnit, EnergyDensityUnit, TemperatureUnit, &
             nValues, VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_I


  subroutine Set_MWV_1_G ( I, Fluid, SpecificOpacity )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

    I % SpecificOpacity  =   SpecificOpacity
    I % Fluid            =>  Fluid

  end subroutine Set_MWV_1_G


  subroutine Set_MWV_1_S &
               ( I, Fluid, Energy, d3_Energy, SpecificOpacity, iBaseCell )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Energy, &
      d3_Energy
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity
    integer ( KDI ), intent ( in ) :: &
      iBaseCell

    I % iBaseCell        =   iBaseCell
    I % SpecificOpacity  =   SpecificOpacity
    I % Energy           =>  Energy
    I % d3_Energy        =>  d3_Energy
    I % Fluid            =>  Fluid

  end subroutine Set_MWV_1_S


  subroutine Compute ( I, R )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    associate &
      (   F => I % Fluid, &
        iBC => I % iBaseCell )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      select type ( R )
      class is ( PhotonMoments_G_Form )
        call I % ComputeKernel_G &
               ( R % Value ( :, R % TEMPERATURE_PARAMETER ), &
                 F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                 F % Value ( :, F % TEMPERATURE ), &
                 I % Value ( :, I % EMISSIVITY_J ), &
                 I % Value ( :, I % OPACITY_J ), &
                 I % Value ( :, I % OPACITY_H ), &
                 I % Value ( :, I % EQUILIBRIUM_J ) )
      end select !-- R
    case ( 'SPECTRAL' )
      call SetPlanckSpectrum &
             ( I % Energy, &
               F % Value ( iBC, F % TEMPERATURE ), &
               I % Value ( :, I % EQUILIBRIUM_J ) )
      call I % ComputeKernel_S &
             ( I % Value ( :, I % EQUILIBRIUM_J ), &
               F % Value ( iBC, F % BARYON_MASS ), &
               F % Value ( iBC, F % COMOVING_BARYON_DENSITY ), &
               F % Value ( iBC, F % TEMPERATURE ), &
               I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ) )
    end select !-- MomentsType

    end associate !-- F, etc.

  end subroutine Compute


  subroutine ComputeTimeScale ( I, R )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      R

    associate &
      (   F => I % Fluid, &
        iBC => I % iBaseCell )

    select case ( trim ( I % MomentsType ) )
    case ( 'GREY' )
      select type ( R )
      class is ( PhotonMoments_G_Form )
        call I % ComputeTimeScaleKernel_G &
               ( R % Value ( :, R % TEMPERATURE_PARAMETER ), &
                 F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                 F % Value ( :, F % INTERNAL_ENERGY ), &
                 F % Value ( :, F % TEMPERATURE ), &
                 R % Value ( :, R % COMOVING_ENERGY ), &
                 I % Value ( :, I % TIME_SCALE ) )
      end select !-- R
    case ( 'SPECTRAL' )
      select type ( R )
      class is ( PhotonMoments_S_Form )
        call SetPlanckSpectrum &
               ( I % Energy, &
                 F % Value ( iBC, F % TEMPERATURE ), &
                 I % Value ( :, I % EQUILIBRIUM_J ) )
        call I % ComputeTimeScaleKernel_S &
               ( I % Value ( :, I % EQUILIBRIUM_J ), &
                 R % Value ( :, R % COMOVING_ENERGY ), &
                 I % d3_Energy, &
                 F % Value ( iBC, F % BARYON_MASS ), &
                 F % Value ( iBC, F % COMOVING_BARYON_DENSITY ), &
                 F % Value ( iBC, F % INTERNAL_ENERGY ), &
                 I % Value ( :, I % TIME_SCALE ) )
      end select !-- R
    end select !-- MomentsType

    end associate !-- F, etc.

  end subroutine ComputeTimeScale


  subroutine ComputeKernel_G ( I, TP, M, N, T, Xi_J, Chi_J, Chi_H, J_EQ )

    class ( Interactions_MWV_1_Form ), intent ( in ) :: &
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
      a, &
      Kappa

    a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa  =  I % SpecificOpacity

    nValues  =  size ( Xi_J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      J_EQ  ( iV )  =  a  *  T ( iV ) ** 4
      Xi_J  ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  J_EQ ( iV )
      Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV ) 
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel_G


  subroutine ComputeKernel_S ( I, J_EQ, M, N, T, Xi_J, Chi_J, Chi_H )

    class ( Interactions_MWV_1_Form ), intent ( in ) :: &
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
      Kappa

    Kappa    =  I % SpecificOpacity
    nValues  =  size ( Xi_J )

    do iV = 1, nValues
      Xi_J  ( iV )  =  Kappa  *  M  *  N  *  J_EQ ( iV )
      Chi_J ( iV )  =  Kappa  *  M  *  N 
      Chi_H ( iV )  =  Chi_J ( iV )
    end do !-- iV

  end subroutine ComputeKernel_S


  subroutine ComputeTimeScaleKernel_G ( I, TP, M, N, E, T, J, TS )

    class ( Interactions_MWV_1_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      TP, &
      M, &
      N, &
      E, &
      T, &
      J
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      TS

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a, &
      Kappa, &
      J_EQ, &
      Q, &
      SqrtTiny

    a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
    Kappa  =  I % SpecificOpacity

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    nValues  =  size ( J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      J_EQ       =  a  *  T ( iV ) ** 4
      Q          =  abs ( Kappa  *  M ( iV )  *  N ( iV )  &
                          *  ( J_EQ - J ( iV ) ) )
      TS ( iV )  =  E ( iV ) / max ( Q, SqrtTiny )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeTimeScaleKernel_G


  subroutine ComputeTimeScaleKernel_S ( I, J_EQ, J, dV, M, N, E, TS )

    class ( Interactions_MWV_1_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J_EQ, &
      J, &
      dV
    real ( KDR ), intent ( in ) :: &
      M, &
      N, &
      E
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      TS

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Kappa, &
      Q, &
      SqrtTiny

    Kappa  =  I % SpecificOpacity

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    nValues  =  size ( J )

    Q = 0.0_KDR
    do iV = 1, nValues
      Q  =  Q  +  abs ( Kappa * M * N * ( J_EQ ( iV ) - J ( iV ) ) ) * dV ( iV )
    end do !-- iV
    TS  =  E / max ( Q, SqrtTiny )

  end subroutine ComputeTimeScaleKernel_S


end module Interactions_MWV_1__Form
