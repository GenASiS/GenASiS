module Interactions_G_G__Form

  !-- Interactions_Generic_Grey_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_G_G_Form
    real ( KDR ) :: &
      OpacityAbsorption = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      Set => Set_G_G
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private, pass :: &
      ComputeKernel
  end type Interactions_G_G_Form

contains


  subroutine InitializeAllocate_I &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_G_G_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_G_G'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_I


  subroutine Set_G_G ( I, OpacityAbsorption )

    class ( Interactions_G_G_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      OpacityAbsorption

    I % OpacityAbsorption  =  OpacityAbsorption

  end subroutine Set_G_G


  subroutine Compute ( I, R, F )

    class ( Interactions_G_G_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( inout ) :: &
      R
    class ( Fluid_P_Template ), intent ( in ) :: &
      F

    call I % ComputeKernel &
           ( F % Value ( :, F % TEMPERATURE ), &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ), &
             I % Value ( :, I % EQUILIBRIUM_J ) )

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_G_G_Form ), intent ( inout ) :: &
      I

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( I, T, Xi_J, Chi_J, Chi_H, J_Eq )

    class ( Interactions_G_G_Form ), intent ( in ) :: &
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

  end subroutine ComputeKernel


end module Interactions_G_G__Form
