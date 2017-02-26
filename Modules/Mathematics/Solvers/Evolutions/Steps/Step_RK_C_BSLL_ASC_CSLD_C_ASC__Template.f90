!-- Step_RK_C_BSLL_ASC_CSLD_C_ASC is a template for a RungeKutta time step of
!   one conserved current on a bundle and another on a base space.

module Step_RK_C_BSLL_ASC_CSLD_C_ASC__Template

  !-- Step_RungeKutta_Current_AtlasSingleChart_1D_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Operations
  use Fields
  use Step_RK_C_ASC__Template
  use Step_RK_C_ASC_1D__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_1D_Template ), abstract :: &
    Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template
      type ( ApplyDivergence_C_Pointer ) :: &
        ApplyDivergenceSection
      type ( ApplySources_C_Pointer ) :: &
        ApplySourcesSection
      type ( ApplyRelaxation_C_Pointer ) :: &
        ApplyRelaxationSection
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC
  end type Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template

contains


  subroutine InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC &
               ( S, NameSuffix, A, B, C )

    class ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_BSLL_ASC_CSLD_C_ASC' 

    call S % InitializeTemplate_C_ASC ( NameSuffix, A, B, C )

  end subroutine InitializeTemplate_C_BSLL_ASC_CSLD_C_ASC


end module Step_RK_C_BSLL_ASC_CSLD_C_ASC__Template
