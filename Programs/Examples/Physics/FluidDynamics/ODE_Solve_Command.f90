module ODE_Solve_Command

  use GenASIS

  implicit none
  private

  type, public :: ODEForm
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaximumSteps = 500000, &
      nSteps       = 0
    real ( KDR ) :: &
      RequestedAccuracy = 1.0e-8_KDR, &
      SolutionAccuracy
    logical ( KDL ) :: &
      Success = .false.
    character ( LDL ) :: &
      Algorithm
    procedure ( DerivativeEvaluatorInterface ), public, pointer, nopass :: &
      FunctionDerivativeEvaluator => null () 
    class ( * ), public, pointer :: &
      FunctionParameters => null ()
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SolveRK4             !-- Fourth Order Runga-Kutta
    procedure, public, pass :: &
      SolveDP5             !-- Dormand-Prince 5th order RK with adaptive steps
    generic :: &
      Solve => SolveRK4, SolveDP5
    final :: &
      Finalize
  end type ODEForm

  abstract interface 
    subroutine DerivativeEvaluatorInterface ( Parameters, X, Y, dYdX )
      use GenASiS
      implicit none
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        Y
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        dYdX
    end subroutine DerivativeEvaluatorInterface
  end interface

    public :: &
      DerivativeEvaluatorInterface, &
      IntegrateODE

contains

  subroutine Initialize &
               ( ODEF, FunctionParameters, FunctionDerivativeEvaluator, &
                 AccuracyOption, MaximumStepOption, VerbosityOption )
    class ( ODEForm ), intent ( inout ) :: &
      ODEF
    class ( * ), intent ( in ), target :: &
      FunctionParameters
    procedure ( DerivativeEvaluatorInterface ), intent ( in ), pointer :: &
      FunctionDerivativeEvaluator
    real ( KDR ), intent ( in ), optional :: &
      AccuracyOption
    integer ( KDI ), intent ( in ), optional :: &
      MaximumStepOption, &
      VerbosityOption
      
    ODEF % IGNORABILITY = CONSOLE % INFO_3 
    if ( present ( VerbosityOption ) ) ODEF % IGNORABILITY = VerbosityOption 
    
    if ( present ( MaximumStepOption ) ) &
      ODEF % MaximumSteps = MaximumStepOption
    
    if ( present ( AccuracyOption ) ) &
      ODEF % RequestedAccuracy = AccuracyOption
    
    ODEF % FunctionParameters          => FunctionParameters
    ODEF % FunctionDerivativeEvaluator => FunctionDerivativeEvaluator

  end subroutine Initialize


  subroutine SolveRK4 ( ODEF, H, X, Y, dYdX, Y_Out ) 
    class ( ODEForm ), intent ( inout ) :: &
      ODEF 
    real ( KDR ), intent ( in ) :: &
      H, &
      X
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y, &
      dYdX
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Y_Out
    

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      HH, &
      H6, &
      XH
    real ( KDR ), dimension ( size ( Y ) ) :: &
      YT, &
      dYT, &
      dYM

    HH = H / 2
    H6 = H / 6.0_KDR
    XH = X + HH

    !-- First step
    do iV = 1, size ( Y ) 
      YT ( iV ) = Y ( iV ) + HH * dYdX ( iV )
    end do

    !-- Second step
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, XH, YT, dYT )
    do iV = 1, size ( Y )
      YT ( iV ) = Y ( iV ) + HH * dYT ( iV )
    end do

    !-- Third Step
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, XH, YT, dYM )
    do iV = 1, size ( Y )
      YT ( iV ) = Y ( iV ) + HH * dYM ( iV )
      dYM ( iV ) = dYT ( iV ) + dYM ( iV )
    end do

    !-- Fourth Step
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, X + H, YT, dYT )
    !-- Accumulate increments with proper weights
    do iV = 1, size ( Y )
      Y_Out ( iV ) = Y ( iV ) &
                       + H6 * ( dYdX ( iV ) + dYT ( iV ) + 2 * dYM ( iV ) )
    end do

  end subroutine SolveRK4

  
  subroutine SolveDP5 ( ODEF, H, X, Y, dYdX, Y_Out, Y_Error ) 
    class ( ODEForm ), intent ( inout ) :: &
      ODEF 
    real ( KDR ), intent ( in ) :: &
      H, &
      X
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y, &
      dYdX
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Y_Out, &
      Y_Error

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      XH
    real ( KDR ), dimension ( size ( Y ) ) :: &
      YT, &
      K2, &
      K3, &
      K4, &
      K5, &
      K6, &
      dYdX_New

    XH = X + H

    !-- First step
    do iV = 1, size ( Y ) 
      YT ( iV ) = Y ( iV ) + H * 0.2 * dYdX ( iV )
    end do

    !-- Second Step
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, X + 0.2 * H, YT, K2 )
    do iV = 1, size ( Y ) 
      YT ( iV ) = Y ( iV ) &
                   + H &
                       * (   3.0_KDR / 40.0_KDR * dYdX ( iV ) &
                           + 9.0_KDR / 40.0_KDR * K2 ( iV ) )
    end do

    !-- Third Step 
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, X + 0.3 * H, YT, K3 )
    do iV = 1, size ( Y ) 
      YT ( iV ) = Y ( iV ) &
                   + H &
                       * (   44.0_KDR / 45.0_KDR * dYdX ( iV ) &
                           - 56.0_KDR / 15.0_KDR * K2 ( iV ) &
                           + 32.0_KDR /  9.0_KDR * K3 ( iV ) )
    end do

    !-- Fourth Step 
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, X + 0.8 * H, YT, K4 )
    do iV = 1, size ( Y ) 
      YT ( iV ) = Y ( iV ) &
                   + H &
                       * (   19372.0_KDR / 6561.0_KDR * dYdX ( iV ) &
                           - 25360.0_KDR / 2187.0_KDR * K2 ( iV ) &
                           + 64448.0_KDR / 6561.0_KDR * K3 ( iV ) &
                           - 212.0_KDR   / 729.0_KDR  * K4 ( iV ) )
    end do

    !-- Fifth Step 
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, X + 8.0_KDR / 9.0_KDR * H, YT, K5 )
    do iV = 1, size ( Y ) 
      YT ( iV ) = Y ( iV ) &
                   + H &
                       * (   9017.0_KDR  / 3168.0_KDR  * dYdX ( iV ) &
                           - 355.0_KDR   / 33.0_KDR    * K2 ( iV ) &
                           + 46732.0_KDR / 5247.0_KDR  * K3 ( iV ) &
                           + 49.0_KDR    / 176.0_KDR   * K4 ( iV ) &
                           - 5103.0_KDR  / 18656.0_KDr * K5 ( iV ))
    end do

    !-- Sixth Step 
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, XH, YT, K6 )
    !-- Accumulate increments with proper weights
    do iV = 1, size ( Y ) 
      Y_Out ( iV ) = Y ( iV ) &
                       + H &
                         * (   35.0_KDR   / 384.0_KDR  * dYdX ( iV ) &
                             + 500.0_KDR  / 1113.0_KDR * K3 ( iV ) &
                             + 125.0_KDR  / 192.0_KDR  * K4 ( iV ) &
                             - 2187.0_KDR / 6784.0_KDr * K5 ( iV ) &
                             + 11.0_KDR   / 84.0_KDR   * K6 ( iV ) )
    end do

    !-- Estimate error as difference between fourth- and fifth- order methods
    call ODEF % FunctionDerivativeEvaluator &
           ( ODEF % FunctionParameters, XH, Y_Out, dYdX_New )

    do iV = 1, size ( Y ) 
      Y_Error ( iV ) = H &
                         * (   71.0_KDR    / 57600.0_KDR  * dYdX ( iV ) &
                             - 71.0_KDR    / 16695.0_KDR  * K3 ( iV ) &
                             + 71.0_KDR    / 1920.0_KDR   * K4 ( iV ) &
                             - 17253.0_KDR / 339200.0_KDr * K5 ( iV ) &
                             + 22.0_KDR    / 525.0_KDR    * K6 ( iV ) &
                             - 1.0_KDR     / 40.0_KDR     * dYdX_New ( iV ) )
    end do

  end subroutine SolveDP5

  
  subroutine IntegrateODE &
               ( ODEF, Y_Start, X1, X2, Epsilon, H1, MaximumStepOption )

    class ( ODEForm ), intent ( inout ) :: &
      ODEF 
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Y_Start
    real ( KDR ), intent ( in ) :: &
      X1, &
      X2, &
      Epsilon, &
      H1
    integer ( KDI ), intent ( in ), optional :: &
      MaximumStepOption

    integer ( KDI ) :: &
      K
    real ( KDR ) :: &
      H_Min, &
      X, &
      H
    real ( KDR ), dimension ( : ), allocatable :: &
      Y, &
      Y_Out, &
      dYdX

    allocate &
      ( Y ( size ( Y_Start ) ), &
        Y_Out ( size ( Y_Start ) ), &
        dYdX ( size ( Y_Start ) ) )

    H_Min = ( X2 - X1 ) / ODEF % MaximumSteps

      Y = Y_Start
      Y_Out = 0.0_KDR
      dYdX = 0.0_KDR
      X = X1
      H = H_Min !-- Take MaximumStep equally sized steps

      do K = 1, ODEF % MaximumSteps
        call ODEF % FunctionDerivativeEvaluator &
                      ( ODEF % FunctionParameters, X, Y, dYdX )
        call ODEF % Solve ( H, X, Y, dYdX, Y_Out ) 
        X = X + H
        Y = Y_Out
      end do

      Y_Start = Y

  end subroutine IntegrateODE

  subroutine Finalize ( ODEF )
    type ( ODEForm ), intent ( inout ) :: &
      ODEF 

    nullify ( ODEF % FunctionParameters )
    nullify ( ODEF % FunctionDerivativeEvaluator )

  end subroutine

end module ODE_Solve_Command
      
