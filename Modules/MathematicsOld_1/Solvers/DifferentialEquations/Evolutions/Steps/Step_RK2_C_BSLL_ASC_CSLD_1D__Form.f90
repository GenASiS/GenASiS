!-- Step_RK2_C_BSLL_ASC_CSLD_1D implements a second-order RungeKutta time 
!   step of an array of conserved currents on a bundle.

module Step_RK2_C_BSLL_ASC_CSLD_1D__Form

  !-- Step_RungeKuttaSecondOrder_Current_BundleSingleLevelDistributed
  !   _AtlasSingleChart_ChartSingleLevelDistributed_Chart_AtlasSingleChart
  !   _Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Fields
  use EvolutionBasics
  use Step_RK_C_BSLL_ASC_CSLD_1D__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ) :: &
    Step_RK2_C_BSLL_ASC_CSLD_1D_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Step_RK2_C_BSLL_ASC_CSLD_1D_Form

contains


  subroutine Initialize ( S, I, Current_BSLL_ASC_CSLD_1D, NameSuffix )

    class ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form ), intent ( inout ) :: &
      S
    class ( IntegratorHeaderForm ), intent ( in ) :: &
      I
    class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
      intent ( in ) :: &
        Current_BSLL_ASC_CSLD_1D
    character ( * ), intent ( in ) :: &
      NameSuffix

    real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
      A
    real ( KDR ), dimension ( 2 : 2 ) :: &
      C
    real ( KDR ), dimension ( 1 : 2 ) :: &
      B

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK2_C_BSLL_ASC_CSLD_1D' 

    call Clear ( A )
    A ( 2, 1 ) = 1.0_KDR

    B ( 1 ) = 0.5_KDR
    B ( 2 ) = 0.5_KDR

    C ( 2 ) = 1.0_KDR
    
    call S % InitializeTemplate_C_BSLL_ASC_CSLD_1D &
           ( I, Current_BSLL_ASC_CSLD_1D, NameSuffix, A, B, C )

  end subroutine Initialize


  impure elemental subroutine Finalize ( S )

    type ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form ), intent ( inout ) :: &
      S

    call S % FinalizeTemplate_C_BSLL_ASC_CSLD_1D ( )

  end subroutine Finalize


end module Step_RK2_C_BSLL_ASC_CSLD_1D__Form
