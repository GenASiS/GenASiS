module MySinusoidalDerivativeFunction_Module
  
  use GenASiS
  
  implicit none
  private

  public :: MySinusoidalDerivativeFunction

contains 

  
  subroutine MySinusoidalDerivativeFunction ( Parameters, X, Y, dYdX )
    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      dYdX

    dYdX ( 1 ) = Y ( 2 )
    dYdX ( 2 ) = - Y ( 1 )
  
  end subroutine MySinusoidalDerivativeFunction


end module MySinusoidalDerivativeFunction_Module

program ODE_Solve_Command_Test

  use GenASiS
  use ODE_Solve_Command
  use MySinusoidalDerivativeFunction_Module

  implicit none

  interface 
    subroutine SinusoidalDerivativeInterface ( Parameters, X, Y, dYdX )
      use GenASiS
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        Y
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        dYdX
    end subroutine SinusoidalDerivativeInterface
  end interface

  type ( ODEForm ), allocatable :: &
    ODEF
  procedure ( SinusoidalDerivativeInterface ), pointer :: &
    Derivative => null ( )
  real ( KDR ), dimension ( 2 ) :: &
    Y_Start, &
    dYdX
  real ( KDR ) :: &
    X1, &
    X2, &
    epsilon, &
    H1, &
    Parameters

  Parameters = 1.0_KDR

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'ODE_Test' )

  Derivative => MySinusoidalDerivativeFunction
  allocate ( ODEF )
  call ODEF % Initialize ( Parameters, Derivative )

  Y_Start ( 1 ) = 2.0_KDR
  Y_Start ( 2 ) = 0.0_KDR

  X1 = 0.0_KDR
  X2 = 1.0_KDR * CONSTANT % PI

  dYdX = 0.0_KDR

  call ODEF % FunctionDerivativeEvaluator &
                ( ODEF % FunctionParameters, X1, Y_Start, dYdX )

  Epsilon = ODEF % RequestedAccuracy

  H1 = ( X2 - X1 ) / ODEF % MaximumSteps

  call ODEF % IntegrateODE ( Y_Start, X1, X2, Epsilon, H1 )

  call Show ( Y_Start ( 1 ), 'Computed X(pi)' )
  call Show ( Y_Start ( 2 ), 'Computed V(pi)' )

  call Show ( 2 * cos ( X2 ), 'Analytic X(pi)' )
  call Show ( -2 * sin ( X2 ), 'Analytic V(pi)')

  call Show ( abs ( ( Y_Start ( 1 ) - 2 * cos ( X2 ) ) / ( 2 * cos ( X2 ) ) ), 'Absolute X Error' )
  call Show ( abs ( ( Y_Start ( 2 ) + 2 * sin ( X2 ) )/ ( 2 * sin ( X2 ) ) ), 'Absolute V Error' )

  deallocate ( PROGRAM_HEADER )

end program ODE_Solve_Command_Test
