module Interactions_C__Form

  !-- Interactions_Constant_Form

  use Basics
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_C_Form
    real ( KDR ) :: &
      EquilibriumDensity = 0.0_KDR, &
      EffectiveOpacity   = 0.0_KDR, &
      TransportOpacity   = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_C
    procedure, public, pass :: &
      Compute
    generic, public :: &
      Initialize => InitializeAllocate_C
    final :: &
      Finalize
  end type Interactions_C_Form

    private :: &
      ComputeKernel

contains


  subroutine InitializeAllocate_C &
               ( I, EquilibriumDensity, EffectiveOpacity, TransportOpacity, &
                 nValues, NameOption, ClearOption, UnitOption )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      EquilibriumDensity, &
      EffectiveOpacity, &
      TransportOpacity
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_C'

    call I % InitializeTemplate ( nValues, NameOption, ClearOption, UnitOption )

    I % EquilibriumDensity  =  EquilibriumDensity
    I % EffectiveOpacity    =  EffectiveOpacity
    I % TransportOpacity    =  TransportOpacity

  end subroutine InitializeAllocate_C


  subroutine Compute ( I )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I

    call ComputeKernel &
           ( I % EquilibriumDensity, I % EffectiveOpacity, &
             I % TransportOpacity, I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
             I % Value ( :, I % EFFECTIVE_OPACITY ), &
             I % Value ( :, I % TRANSPORT_OPACITY ) )

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_C_Form ), intent ( inout ) :: &
      I

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeKernel ( ED, EO, TO, EDV, EOV, TOV )

    real ( KDR ), intent ( in ) :: &
      ED, &
      EO, &
      TO
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      EDV, &
      EOV, &
      TOV

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues  =  size ( EDV )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      EDV ( iV )  =  EDV ( iV ) ! ED
      EOV ( iV )  =  EOV ( iV ) ! EO
      TOV ( iV )  =  TOV ( iV ) ! TO
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeKernel


end module Interactions_C__Form
