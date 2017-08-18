!-- Step_RK_C_BSLL_ASC_CSLD_C_ASC implements a second-order RungeKutta time 
!   step of one conserved current on a bundle and another on its base space.

module Step_RK2_C_BSLL_ASC_CSLD_C_ASC__Form

  !-- Step_RungeKuttaSecondOrder_Current_BundleSingleLevelDistributed
  !   _AtlasSingleChart_ChartSingleLevelDistributed_Chart_AtlasSingleChart
  !   _Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Fields
  use Step_RK_C_BSLL_ASC_CSLD_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ) :: &
    Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form

contains


  subroutine Initialize ( S, Current_BSLL_ASC_CSLD, Current_ASC, NameSuffix )

    class ( Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form ), intent ( inout ) :: &
      S
    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ) :: &
      Current_BSLL_ASC_CSLD
    class ( Current_ASC_Template ), intent ( in ) :: &
      Current_ASC
    character ( * ), intent ( in ) :: &
      NameSuffix

    real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
      A
    real ( KDR ), dimension ( 2 : 2 ) :: &
      C
    real ( KDR ), dimension ( 1 : 2 ) :: &
      B

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK2_C_BSLL_ASC_CSLD_C_ASC' 

    call Clear ( A )
    A ( 2, 1 ) = 1.0_KDR

    B ( 1 ) = 0.5_KDR
    B ( 2 ) = 0.5_KDR

    C ( 2 ) = 1.0_KDR
    
    call S % InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC &
           ( Current_BSLL_ASC_CSLD, Current_ASC, NameSuffix, A, B, C )

  end subroutine Initialize


  impure elemental subroutine Finalize ( S )

    type ( Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form ), intent ( inout ) :: &
      S

    call S % FinalizeTemplate_C_BSLL_ASC_CSLD_C_ASC ( )

  end subroutine Finalize


end module Step_RK2_C_BSLL_ASC_CSLD_C_ASC__Form
