!-- Step_RK_C_BSLL_ASC_CSLD_1D is a template for a RungeKutta time step of
!   an array of conserved currents on a bundle.

module Step_RK_C_BSLL_ASC_CSLD_1D__Template

  !-- Step_RungeKutta_Current_BundleSingleLevelDistributed_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Template

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Manifolds
  use Operations
  use Fields
  use EvolutionBasics
  use Increments
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ), abstract :: &
    Step_RK_C_BSLL_ASC_CSLD_1D_Template
      integer ( KDI ) :: &
        iTimerSections = 0, &
        iTimerStoreSections = 0, &
        iTimerFibers = 0, &
        iTimerLoadSections = 0, &
        nFibers   = 0, &
        nSections = 0, &
        nCurrents = 0
      type ( Real_3D_2D_Form ), dimension ( :, : ), allocatable :: &
        BoundaryFluence_CSL_S
      class ( Chart_SLL_Form ), pointer :: &
        Chart_F => null ( )
      class ( Bundle_SLL_ASC_CSLD_Form ), pointer :: &
        Bundle_SLL_ASC_CSLD => null ( )
      class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
        pointer :: &
          Current_BSLL_ASC_CSLD_1D => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate_C_BSLL_ASC_CSLD_1D
  end type Step_RK_C_BSLL_ASC_CSLD_1D_Template

contains


  subroutine InitializeTemplate_C_BSLL_ASC_CSLD_1D &
               ( S, I, Current_BSLL_ASC_CSLD_1D, NameSuffix, A, B, C )

    class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), intent ( inout ) :: &
      S
    class ( IntegratorHeaderForm ), intent ( in ) :: &
      I
    class ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
      intent ( in ), target :: &
        Current_BSLL_ASC_CSLD_1D
    character ( * ), intent ( in ) :: &
      NameSuffix
    real ( KDR ), dimension ( 2 : , : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( 2 : ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iS, &  !-- iSection
      iC     !-- iCurrent

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK_C_BSLL_ASC_CSLD_C_ASC_1D' 

    call S % InitializeTemplate ( I, NameSuffix, A, B, C )

    S % Current_BSLL_ASC_CSLD_1D  =>  Current_BSLL_ASC_CSLD_1D
    S % nCurrents = size ( Current_BSLL_ASC_CSLD_1D )
    call Show ( S % nCurrents, 'nCurrents', S % IGNORABILITY )

    S % Bundle_SLL_ASC_CSLD  &
      =>  Current_BSLL_ASC_CSLD_1D ( 1 ) % Element % Bundle_SLL_ASC_CSLD

    associate ( B => S % Bundle_SLL_ASC_CSLD )
    S % Chart      =>  B % Base_CSLD
    S % nFibers    =   B % nFibers
    S % nSections  =   B % nSections
    S % Chart_F    =>  B % Fiber_CSLL
    end associate !-- B
    call Show ( S % nFibers,   'nFibers',   S % IGNORABILITY )
    call Show ( S % nSections, 'nSections', S % IGNORABILITY )

    S % Allocated = .false.

  end subroutine InitializeTemplate_C_BSLL_ASC_CSLD_1D


end module Step_RK_C_BSLL_ASC_CSLD_1D__Template
