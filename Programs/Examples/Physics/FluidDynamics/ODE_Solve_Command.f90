module ODE_Solve_Command

  use GenASIS

  implicit none
  private

  type, public :: ODEForm
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaximumSteps = 500000
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
    procedure, public, pass :: &
      IntegrateODE
    procedure, public, pass :: &
      VariableStepSolve
    procedure, public, pass :: &
      VariableStepDoPr
    final :: &
      Finalize
  end type ODEForm

    real ( KDR ), PARAMETER :: &
      P_Grow    = - 0.2_KDR, &
      P_Shrink  = - 0.25_KDR, &
      F_Cor     = 1.0_KDR / 15.0_KDR, &
      Safety    = 0.9_KDR, &
      ErrorCon  = ( 4.0_KDR / Safety ) ** ( 1 / P_Grow ), &
      beta      = 0.4_KDR / 5.0_KDR, &
      alpha     = - ( 0.2_KDR - 0.75_KDR * beta ), &
      MinScale = 0.2_KDR, &
      MaxScale = 10.0_KDR



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
      DerivativeEvaluatorInterface

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


  subroutine Finalize ( ODEF )
    type ( ODEForm ), intent ( inout ) :: &
      ODEF 

    nullify ( ODEF % FunctionParameters )
    nullify ( ODEF % FunctionDerivativeEvaluator )

  end subroutine


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
               ( ODEF, Y_Start, X1, X2, Epsilon, H1, RK4Option )

    class ( ODEForm ), intent ( inout ) :: &
      ODEF 
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Y_Start
    real ( KDR ), intent ( in ) :: &
      X1, &
      X2, &
      Epsilon, &
      H1
    logical ( KDL ), intent ( in ), optional :: &
      RK4Option

    logical ( KDL ) :: &
      RK4, &
      Reject
    integer ( KDI ) :: &
      nS, &             !-- number of steps
      nOk, &            !-- number of successful steps taken
      nBad              !-- number of unsuccessful steps taken
    real ( KDR ) :: &
      H_Min, &
      X, &
      H, &
      H_Out, &
      H_New, &
      Error
    real ( KDR ), dimension ( : ), allocatable :: &
      Y, &
      Y_Scale, &
      dYdX

    allocate &
      ( Y ( size ( Y_Start ) ), &
        Y_Scale ( size ( Y_Start ) ), &
        dYdX ( size ( Y_Start ) ) )

    X = X1
    H = H1
    Y = Y_Start

    H_Out = 0.0_KDR
    H_New = 0.0_KDR
    Error = 1.0e-4_KDR

    nOk  = 0
    nBad = 0

    RK4 = .true. 
    if ( present ( RK4Option ) ) &
      RK4 = RK4Option

    Reject = .true.

    do nS = 1, ODEF % MaximumSteps
      call ODEF % FunctionDerivativeEvaluator &
                      ( ODEF % FunctionParameters, X, Y, dYdX )
      if ( X + H > X2 ) &
          H = X2 - X
      if ( RK4 ) then
        Y_Scale = abs ( Y + abs ( H * dYdX ) ) + sqrt ( tiny ( 0.0_KDR ) )
        call ODEF % VariableStepSolve &
                    ( Y, dYdX, X, H, Epsilon, Y_Scale, H_Out, H_New )
      else
        call ODEF % VariableStepDoPr &
                      ( Y, dYdX, X, Error, Reject, H, &
                        Epsilon, H_Out, H_New )
      end if

      if ( H_Out == H ) then 
        nOk = nOk + 1
      else
        nBad = nBad + 1
      end if
      
      if ( X >= X2 ) then
        Y_Start = Y
        call Show ( 'Solution obtained in')
        call Show ( nS, 'total steps' )
        call Show ( nOk, 'good steps' )
        call Show ( nBad, 'altered steps' )
        return
      end if

      H = H_New
    end do

    call Show ( 'Solution obtained in maximum number of steps allowed' )

  end subroutine IntegrateODE


  subroutine VariableStepSolve &
               ( ODEF, Y, dYdX, X, H_In, Epsilon_In, Y_Scale, H_Out, H_New )

    class ( ODEForm ), intent ( inout ) :: &
      ODEF 
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Y, &
      dYdX
    real ( KDR ), intent ( inout ) :: &
      X
    real ( KDR ), intent ( in ) :: &
      H_In, &
      Epsilon_In
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y_Scale
    real ( KDR ), intent ( out ) :: &
      H_Out, &
      H_New


    real ( KDR ) :: &
      H_Min, &
      X_Save, &
      H, &
      HH
    real ( KDR ), dimension ( size ( Y ) ) :: &
      Y_Save, &
      Y_Temp, &
      dYdX_Save, &
      Epsilon, &
      ErrorMax 

    X_Save    = X
    Y_Save    = Y
    dYdX_Save = dYdX

    H = H_In
    HH = H / 2

    call ODEF % Solve ( HH, X_Save, Y_Save, dYdX_Save, Y_Temp ) 

    X = X_Save + HH

    call ODEF % FunctionDerivativeEvaluator &
                  ( ODEF % FunctionParameters, X, Y_Temp, dYdX )

    call ODEF % Solve ( HH, X, Y_Temp, dYdX, Y )

    X = X_Save + H
    
    call ODEF % Solve ( H, X_Save, Y_Save, dYdX_Save, Y_Temp ) 

    ErrorMax = 0.0_KDR

    Y_Temp = Y - Y_Temp

    ErrorMax = Max ( ErrorMax, abs ( Y_Temp / Y_Scale ) ) 
    Epsilon = Epsilon_In * Y_Scale

    ErrorMax = ErrorMax / Epsilon

    if ( maxval ( ErrorMax ) > 1.0_KDR ) then
      H = Safety * H * ( maxval ( ErrorMax ) ** P_Shrink )
      call ODEF % VariableStepSolve &
                    ( Y_Save, dYdX_Save, X_Save, H, &
                      Epsilon_In, Y_Scale, H_Out, H_New )
    else
      H_Out = H
      if ( maxval ( ErrorMax ) > ErrorCon ) then
        H_New = Safety * H * ( maxval ( ErrorMax ) ** P_Grow )
      else
        H_New = 4.0_KDR * H
      end if
      Y = Y + Y_Temp * F_Cor
    end if

  end subroutine VariableStepSolve
  

  subroutine VariableStepDoPr & 
               ( ODEF, Y, dYdX, X, Error_Old, Reject, &
                 H_In, Epsilon_In, H_Out, H_New )

    class ( ODEForm ), intent ( inout ) :: &
      ODEF 
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Y, &
      dYdX
    real ( KDR ), intent ( inout ) :: &
      X, &
      Error_Old
    logical ( KDL ), intent ( inout ) :: &
      Reject
    real ( KDR ), intent ( in ) :: &
      H_In, &
      Epsilon_In
    real ( KDR ), intent ( out ) :: &
      H_Out, &
      H_New

    real ( KDR ) :: &
      X_Save, &
      H, &
      Error
    real ( KDR ), dimension ( size ( Y ) ) :: &
      Y_Save, &
      dYdX_Save, &
      Y_Temp, &
      Y_Error

    H = H_In

    call ODEF % Solve ( H, X, Y, dYdX, Y_Temp, Y_Error ) 

    Error &
      = Epsilon_In &
          * ( 1.0_KDR &
              + max ( maxval ( abs ( Y ) ), maxval ( abs ( Y_Temp ) ) ) ) 

    Error = sqrt ( maxval ( Y_Error ) / Error )

    call EvaluateError ( Error, Error_Old, H, H_New, Reject )

    if ( Reject ) then
      call ODEF % VariableStepDoPr &
                    ( Y, dYdX, X, Error, Reject, &
                      H, Epsilon_In, H_Out, H_New )
    else
     H_Out = H
     Y     = Y_Temp
     X     = X + H
    end if

  end subroutine VariableStepDoPr


  subroutine EvaluateError &
               ( Error, Error_Old, H, H_New, Reject )
    real ( KDR ), intent ( inout ) :: &
      Error, &
      Error_Old, &
      H, &
      H_New
    logical ( KDL ), intent ( inout ) :: &
      Reject

    real ( KDR ) :: &
      Scale

    if ( Error <= 1.0_KDR ) then
      if ( Error == 0.0_KDR ) then
        Scale = MaxScale
      else
        Scale = Safety * Error ** alpha * Error_Old ** beta
        if ( Scale < MinScale ) &
          Scale = MinScale
        if ( Scale > MaxScale ) &
          Scale = MaxScale
        if ( Reject ) then
          H_New = H * min ( Scale, 1.0_KDR )
        else
          H_New = H * Scale
          Error_Old = max ( Error, 1.0e-4_KDR )
        end if
      end if
      Reject = .false.
    else
      Scale = max ( Safety * Error ** alpha, MinScale )
      H = H * Scale
      Reject = .true.
    end if

  end subroutine EvaluateError

end module ODE_Solve_Command
      
