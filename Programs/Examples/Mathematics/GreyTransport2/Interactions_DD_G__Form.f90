module Interactions_DD_G__Form

  !-- Interactions_DynamicDiffusion_Grey__Form, Just et al. 2015
  
  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_DD_G_Form
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR, &
      D, &
      Delta, &
      X0, &
      Y0
    real ( KDR ), pointer :: &
      Time
    class ( GeometryFlatForm ), pointer :: &
      Geometry => null ( ) 
    class ( RadiationMomentsForm ), pointer :: &
      Radiation => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_DD_G
    generic, public :: &
      Initialize => InitializeAllocate_DD_G
    procedure, private, pass :: &
      Set_DD_G
    generic, public :: &
      Set => Set_DD_G
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private :: &
      ComputeKernel
  end type Interactions_DD_G_Form

contains


  subroutine InitializeAllocate_DD_G &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_DD_G_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_DD_1_G'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_DD_G


  subroutine Set_DD_G &
               ( I, Radiation, Geometry, Time, &
                 D, Delta, X0, Y0 )

    class ( Interactions_DD_G_Form ), intent ( inout ) :: &
      I
    class ( RadiationMomentsForm ), intent ( in ), target :: &
      Radiation
    class ( GeometryFlatForm ), intent ( in ), target :: &
      Geometry
    real ( KDR ), intent ( in ), target :: &
      Time
    real ( KDR ), intent ( in ) :: &
      D, &
      Delta, &
      X0, &
      Y0

    I % Radiation        =>  Radiation
    I % Geometry         =>  Geometry
    I % Time             =>  Time
    

    I % D     = D
    I % Delta = Delta
    I % X0    = X0 
    I % Y0    = Y0

  end subroutine Set_DD_G


  subroutine Compute ( I, Current )

    class ( Interactions_DD_G_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    associate &
      ( R => I % Radiation, &
        G => I % Geometry )

    call I % ComputeKernel &
           ( I % Time, &
             I % D, &
             I % Delta, &
             I % X0, &
             I % Y0, &
             G % Value ( :, G % CENTER ( 1 ) ), &
             G % Value ( :, G % CENTER ( 2 ) ), &
             R % Value ( :, R % FLUID_VELOCITY_U ( 1 ) ), &
             R % Value ( :, R % FLUID_VELOCITY_U ( 2 ) ), &
             R % Value ( :, R % FLUID_VELOCITY_U ( 3 ) ), &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ) )

    end associate !-- R, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_DD_G_Form ), intent ( inout ) :: &
      I

    nullify ( I % Radiation )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel &
               ( I, Time, D, Delta, X0, Y0, &
                 X, Y, V_1, V_2, V_3, Xi_J, Chi_J, Chi_H )

    class ( Interactions_DD_G_Form ), intent ( in ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      Time, &
      D, &
      Delta, &
      X0, &
      Y0
    real ( KDR ), dimension ( : ), intent ( in ) ::  &
      X, &
      Y, &
      V_1, &
      V_2, &
      V_3
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      R_sq

    nValues  =  size ( Xi_J )
    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      R_sq          =   ( X ( iV ) - V_1 ( iV ) * Time - X0 ) ** 2 &
                      + ( Y ( iV ) - V_2 ( iV ) * Time - Y0 ) ** 2
      Xi_J  ( iV )  =  0.0_KDR
      Chi_J ( iV )  =  0.0_KDR 
      Chi_H ( iV )  =  1 / ( 3.0_KDR * D ) * exp ( - R_sq / delta ** 2 )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_DD_G__Form
