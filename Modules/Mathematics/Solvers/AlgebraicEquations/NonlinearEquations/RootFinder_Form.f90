module RootFinder_Form

  use Basics
   
  implicit none
  private
  
  type, public :: RootFinderForm 
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaxIterations = 20, &
      nIterations      = 0
    real ( KDR ) :: &
      RequestedAccuracy = 1.0e-8_KDR, &
      SolutionAccuracy
    logical ( KDL ) :: &
      Success = .false.
    procedure ( EZ ), public, pointer, nopass :: &
      EvaluateZero => null ( ), &
      EvaluateDerivative => null ( )
    class ( * ), private, pointer :: &
      Parameters => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SolveSecant
    procedure, private, pass :: &
      Solve_B_NR
    procedure, private, pass :: &
      SolveBrent
    procedure, private, pass :: &
      SolveNewtonRaphson
    generic :: &
      Solve => SolveSecant, Solve_B_NR
    final :: &
      Finalize
  end type RootFinderForm
  
  abstract interface
    subroutine EZ ( Parameters, Input, Result )
      use Basics
      implicit none
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        Input
      real ( KDR ), intent ( out ) :: &
        Result
    end subroutine EZ
  end interface
  
contains


  subroutine Initialize &
               ( RF, Parameters, AccuracyOption, MaxIterationsOption, &
                 VerbosityOption )

    class ( RootFinderForm ), intent ( inout ) :: &
      RF
    class ( * ), intent ( in ), target :: &
      Parameters
    real ( KDR ), intent ( in ), optional :: &
      AccuracyOption
    integer ( KDI ), intent ( in ), optional :: &
      MaxIterationsOption, &
      VerbosityOption
      
    RF % IGNORABILITY = CONSOLE % INFO_3 
    if ( present ( VerbosityOption ) ) RF % IGNORABILITY = VerbosityOption 
    
    if ( present ( MaxIterationsOption ) ) &
      RF % MaxIterations = MaxIterationsOption
    
    if ( present ( AccuracyOption ) ) &
      RF % RequestedAccuracy = AccuracyOption
    
    RF % Parameters => Parameters
  
  end subroutine Initialize
  
  
  subroutine SolveSecant ( RF, Guess_1, Guess_2, Root )
  
    class ( RootFinderForm ), intent ( inout ) :: &
      RF
    real ( KDR ), intent ( in ) :: &
      Guess_1, &
      Guess_2
    real ( KDR ), intent ( out ) :: &
      Root
      
    integer ( KDI ) :: &
      iIteration
    real ( KDR ) :: &
      X, &
      X_0, X_1, &
      Y_0, Y_1
      
    RF % Success = .false.
    RF % nIterations = 0
      
    X_0 = Guess_1
    X_1 = Guess_2
    
    call RF % EvaluateZero ( RF % Parameters, X_0, Y_0 )
    call RF % EvaluateZero ( RF % Parameters, X_1, Y_1 )
    
    do iIteration = 1, RF % MaxIterations
      
      RF % nIterations = iIteration
      
      RF % SolutionAccuracy = abs ( Y_1 )
      Root = X_1
      
      if ( abs ( ( X_1 - X_0 ) /  X_0 ) <= RF % RequestedAccuracy &
           .or. abs ( Y_1 ) <= RF % RequestedAccuracy ) &
      then
        RF % Success = .true. 
        exit
      end if
      
      if ( Y_1 == Y_0 ) &
        exit

      X = X_1 - Y_1 * ( X_1 - X_0 ) / ( Y_1 - Y_0 )
      
      X_0 = X_1
      Y_0 = Y_1
      
      X_1 = X
      call RF % EvaluateZero ( RF % Parameters, X_1, Y_1 )
      
    end do
    
  end subroutine SolveSecant
  

  subroutine Solve_B_NR ( RF, Interval, Root )

    class ( RootFinderForm ), intent ( inout ) :: &
      RF
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      Interval
    real ( KDR ), intent ( out ) :: &
      Root
    
    if ( associated ( RF % EvaluateDerivative ) ) then
      call RF % SolveNewtonRaphson ( Interval, Root )
    else
      call RF % SolveBrent ( Interval, Root )
    end if
  end subroutine Solve_B_NR


  subroutine SolveBrent ( RF, Interval, Root )
  
    class ( RootFinderForm ), intent ( inout ) :: &
      RF
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      Interval
    real ( KDR ), intent ( out ) :: &
      Root
      
    integer ( KDI ) :: &
      iIteration
    real ( KDR ) :: &
      a,b,c,d,e,fa,fb,fc,p,q,r,s,Accuracy1,xm
    
    !-- Van Wijngaarden-Dekker-Brent Method, based on Numerical Recipes in
    !   Fortran (1992) routine "zbrent"    
    
    RF % Success = .false.
    
    a = Interval ( 1 )
    b = Interval ( 2 )
    
    call RF % EvaluateZero ( RF % Parameters, a, fa )
    call RF % EvaluateZero ( RF % Parameters, b, fb )
    
    if ( ( fa > 0.0_KDR .and. fb > 0.0_KDR ) &
         .or. ( fa < 0.0_KDR .and. fb < 0.0_KDR ) ) &
    then
      call Show &
             ( 'Invalid interval was given as arguments', RF % IGNORABILITY )
      call Show ( 'RootFinderForm', 'Class', RF % IGNORABILITY )
      call Show ( 'Solve', 'Method', RF % IGNORABILITY )
      call Show ( Interval, 'Interval', RF % IGNORABILITY )
      RF % Success = .false.
      return
    end if
    
    c = b
    fc = fb
    do iIteration = 1, RF % MaxIterations 
      RF % nIterations = iIteration
      if ((fb > 0.0_KDR .and. fc > 0.0_KDR ) &
            .or. (fb < 0.0_KDR .and.fc < 0.0_KDR)) &
      then
        c=a
        fc=fa
        d=b-a
        e=d
      end if
      if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
      end if
      Accuracy1=2.0_KDR*epsilon(1.0_KDR)*abs(b)+0.5_KDR*RF % RequestedAccuracy
      xm=0.5_KDR*(c-b)
      if (abs(xm) <= Accuracy1 .or. fb == 0.0_KDR) then
        Root=b
        RF % Success = .true.
        RF % SolutionAccuracy = abs ( xm )
        return
      end if
      if (abs(e) >= Accuracy1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
          p=2.0_KDR*xm*s
          q=1.0_KDR-s
        else
          q=fa/fc
          r=fb/fc
          p=s*(2.0_KDR*xm*q*(q-r)-(b-a)*(r-1.0_KDR))
          q=(q-1.0_KDR)*(r-1.0_KDR)*(s-1.0_KDR)
        end if
        if (p > 0.0_KDR) q=-q
        p=abs(p)
        if (2.0_KDR*p  <  min(3.0_KDR*xm*q-abs(Accuracy1*q),abs(e*q))) then
          e=d
          d=p/q
        else
          d=xm
          e=d
        end if
      else
        d=xm
        e=d
      end if
      a=b
      fa=fb
      b=b+merge(d,sign(Accuracy1,xm), abs(d) > Accuracy1 )
      call RF % EvaluateZero ( RF % Parameters, b, fb )
    end do
    
    RF % SolutionAccuracy = abs ( fb )
    Root=b
    
    RF % Success = .false.

  end subroutine SolveBrent
  
  
  subroutine SolveNewtonRaphson ( RF, Interval, Root )
    
    !-- Newton-Raphson Method, based on Numerical Recipes in
    !   Fortran (1992) routine "rtsafe" 
  
    class ( RootFinderForm ), intent ( inout ) :: &
      RF
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      Interval
    real ( KDR ), intent ( out ) :: &
      Root
      
    integer ( KDI ) :: &
      iIteration
    real ( KDR ) :: &
      x1, x2, xacc, &
      df, dx, dxold, f, fh, fl, temp, xh, xl

    xacc = RF % RequestedAccuracy ! epsilon ( 1.0_KDR )  * 1.0e1_KDR
    
    x1 = Interval(1)
    x2 = Interval(2)

    !call Function(Parameters, x1, fl)
    !call Derivative(Parameters, x1, df)
    !call Function(Parameters, x2, fh)
    !call Derivative(Parameters, x2, df)
    
    call RF % EvaluateZero ( RF % Parameters, x1, fl )
    call RF % EvaluateDerivative ( RF % Parameters, x1, df ) 
    call RF % EvaluateZero ( RF % Parameters, x2, fh )
    call RF % EvaluateDerivative ( RF % Parameters, x2, df ) 
    
    if ( ( fl > 0.0_KDR .and. fh > 0.0_KDR ) &
         .or. ( fl < 0.0_KDR .and. fh < 0.0_KDR ) ) then
      call Show &
             ( 'Invalid interval was given as arguments', RF % IGNORABILITY )
      call Show ( 'RootFinderForm', 'Class', RF % IGNORABILITY )
      call Show ( 'Solve', 'Method', RF % IGNORABILITY )
      call Show ( Interval, 'Interval', RF % IGNORABILITY )
      RF % Success = .false.
      return
    end if
   
    !if(present(SuccessOption)) SuccessOption = .true.
    
    if ( fl == 0.0_KDR ) then
      Root=x1
      return
    else if ( fh == 0.0_KDR ) then
      Root=x2
      return
    else if ( fl < 0.0_KDR ) then
      xl=x1
      xh=x2
    else
      xh=x1
      xl=x2
    endif
    Root  = 0.5_KDR * ( x1 + x2 )
    dxold = abs ( x2 - x1 )
    dx    = dxold
    
    !call Function(Parameters, Root, f)
    !call Derivative(Parameters, Root, df)

    call RF % EvaluateZero ( RF % Parameters, Root, f )
    call RF % EvaluateDerivative ( RF % Parameters, Root, df ) 
    
    do iIteration = 1, RF % MaxIterations
      if ( ( ( Root-xh ) * df - f ) * ( ( Root - xl ) * df - f ) > 0.0_KDR &
           .or. abs ( 2.0_KDR * f ) > abs ( dxold * df ) ) then
        dxold = dx
        dx    = 0.5_KDR * ( xh - xl )
        Root  = xl + dx
        if ( xl == Root ) return
      else
        dxold = dx
        dx    = f/df
        temp  = Root
        Root  = Root - dx
        if ( temp == Root )return
      endif
      if ( abs ( dx ) < xacc ) return
      call RF % EvaluateZero ( RF % Parameters, Root, f )
      call RF % EvaluateDerivative ( RF % Parameters, Root, df ) 
      if ( f < 0.0_KDR ) then
        xl = Root
      else
        xh = Root
      endif
    end do
    
    call Show ( 'FindRoot exceeded maximum iterations', RF % IGNORABILITY )
    call Show ( 'FindRoot could not find root to the specified accuracy', &
                RF % IGNORABILITY )

     RF % Success = .false.

  end subroutine SolveNewtonRaphson

  
  Subroutine Finalize ( RF )
  
    type ( RootFinderForm ), intent ( inout ) :: &
      RF
    
    nullify ( RF % Parameters )
    nullify ( RF % EvaluateDerivative )
    nullify ( RF % EvaluateZero )
  
  end subroutine Finalize
  

end module RootFinder_Form
