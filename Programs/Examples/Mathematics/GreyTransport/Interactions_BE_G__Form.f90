module Interactions_BE_G__Form

  !-- Interactions_BoseEinstein_Grey__Form
  
  use Basics
  use Fluid_P__Template
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_BE_G_Form
    real ( KDR ) :: &
      SpecificOpacity
    class ( Fluid_P_Template ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_BE_G
    generic, public :: &
      Initialize => InitializeAllocate_BE_G
    procedure, public, pass :: &
      Set
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Interactions_BE_G_Form

    private :: &
      ComputeKernel

contains


  subroutine InitializeAllocate_BE_G &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_BE_G_Form ), intent ( inout ) :: &
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
      I % Type = 'an Interactions_BE_G'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, nValues, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_BE_G


  subroutine Set ( I, Fluid, SpecificOpacity )

    class ( Interactions_BE_G_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Template ), intent ( in ), target :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

    I % SpecificOpacity  =  SpecificOpacity
    I % Fluid  =>  Fluid

  end subroutine Set


  subroutine Compute ( I )

    class ( Interactions_BE_G_Form ), intent ( inout ) :: &
      I

    associate ( F => I % Fluid )
    call ComputeKernel &
           ( F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % COMOVING_DENSITY ), &
             F % Value ( :, F % TEMPERATURE ), &
             I % SpecificOpacity, CONSTANT % RADIATION, &
             I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
             I % Value ( :, I % EFFECTIVE_OPACITY ), &
             I % Value ( :, I % TRANSPORT_OPACITY ) )
    end associate !-- F

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_BE_G_Form ), intent ( inout ) :: &
      I

    nullify ( I % Fluid )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( M, N, T, Kappa, a, EDV, EOV, TOV )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T
    real ( KDR ), intent ( in ) :: &
      Kappa, &
      a
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EDV, &
      EOV, &
      TOV

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues  =  size ( EDV )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      EDV ( iV )  =  a  *  T ( iV ) ** 4
      EOV ( iV )  =  Kappa  *  M ( iV )  *  N ( iV ) 
      TOV ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_BE_G__Form
