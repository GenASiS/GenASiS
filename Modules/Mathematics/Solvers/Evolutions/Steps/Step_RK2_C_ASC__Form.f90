!-- Step_RK2_C_ASC implements a second-order RungeKutta time step of a
!   conserved current.

module Step_RK2_C_ASC__Form

  !-- Step_RungeKuttaSecondOrder_Current_AtlasSingleChart_Form

  !-- See Wikipedia "Runge-Kutta methods" for explanation of Butcher 
  !   tableau entries A, B, C

  use Basics
  use Step_RK_C_ASC__Template

  implicit none
  private

  type, public, extends ( Step_RK_C_ASC_Template ) :: Step_RK2_C_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Step_RK2_C_ASC_Form

contains


  subroutine Initialize ( S, NameSuffix )

    class ( Step_RK2_C_ASC_Form ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ) :: &
      NameSuffix

    real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
      A
    real ( KDR ), dimension ( 2 : 2 ) :: &
      C
    real ( KDR ), dimension ( 1 : 2 ) :: &
      B

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK2_C' 

    call Clear ( A )
    A ( 2, 1 ) = 1.0_KDR

    B ( 1 ) = 0.5_KDR
    B ( 2 ) = 0.5_KDR

    C ( 2 ) = 1.0_KDR
    
    call S % InitializeTemplate_C ( NameSuffix, A, B, C )

  end subroutine Initialize


  impure elemental subroutine Finalize ( S )

    type ( Step_RK2_C_ASC_Form ), intent ( inout ) :: &
      S

    call S % FinalizeTemplate_C ( )

  end subroutine Finalize


end module Step_RK2_C_ASC__Form
