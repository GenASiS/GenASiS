module RootFinding_Form

  use Basics
   
  implicit none
  private
  
  type, public :: RootFindingForm 
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaximumIteration = 100, &
      nIterations      = 0
    real ( KDR ) :: &
      RequestedAccuracy = 1.0e-8_KDR, &
      SolutionAccuracy
    logical ( KDL ) :: &
      Success = .false.
    character ( LDL ) :: &
      Algorithm
    procedure ( ZeroFunctionEvaluatorInterface ), private, pointer, nopass :: &
      FunctionEvaluator => null ( ), &
      FunctionDerivativeEvaluator  => null ( )
    class ( * ), private, pointer :: &
      FunctionParameters => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SolveBrent
    procedure, private, pass :: &
      SolveSecant
    generic :: &
      Solve => SolveBrent, SolveSecant
    final :: &
      Finalize
  end type RootFindingForm
  
  abstract interface
    subroutine ZeroFunctionEvaluatorInterface ( Parameters, Input, Result )
      use Basics
      implicit none
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        Input
      real ( KDR ), intent ( out ) :: &
        Result
    end subroutine ZeroFunctionEvaluatorInterface
  end interface
  
  public :: &
    ZeroFunctionEvaluatorInterface
  
contains


  subroutine Initialize &
               ( RF, FunctionParameters, FunctionEvaluator, AccuracyOption, &
                 MaximumIterationOption, VerbosityOption )
    class ( RootFindingForm ), intent ( inout ) :: &
      RF
    class ( * ), intent ( in ), target :: &
      FunctionParameters
    procedure ( ZeroFunctionEvaluatorInterface ), intent ( in ), pointer :: &
      FunctionEvaluator
    real ( KDR ), intent ( in ), optional :: &
      AccuracyOption
    integer ( KDI ), intent ( in ), optional :: &
      MaximumIterationOption, &
      VerbosityOption
      
    RF % IGNORABILITY = CONSOLE % INFO_3 
    if ( present ( VerbosityOption ) ) RF % IGNORABILITY = VerbosityOption 
    
    if ( present ( MaximumIterationOption ) ) &
      RF % MaximumIteration = MaximumIterationOption
    
    if ( present ( AccuracyOption ) ) &
      RF % RequestedAccuracy = AccuracyOption
    
    RF % FunctionParameters => FunctionParameters
    RF % FunctionEvaluator => FunctionEvaluator
  
  end subroutine Initialize
  
  
  subroutine SolveBrent ( RF, Interval, Root )
  
    class ( RootFindingForm ), intent ( inout ) :: &
      RF
    real ( KDR ), dimension ( : ), intent ( in ) :: &
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
    
    call RF % FunctionEvaluator ( RF % FunctionParameters, a, fa )
    call RF % FunctionEvaluator ( RF % FunctionParameters, b, fb )
    
    if ( ( fa > 0.0_KDR .and. fb > 0.0_KDR ) &
         .or. ( fa < 0.0_KDR .and. fb < 0.0_KDR ) ) &
    then
      call Show &
             ( 'Invalid interval was given as arguments', RF % IGNORABILITY )
      call Show ( 'RootFindingForm', 'Class', RF % IGNORABILITY )
      call Show ( 'Solve', 'Method', RF % IGNORABILITY )
      call Show ( Interval, 'Interval', RF % IGNORABILITY )
      RF % Success = .false.
      return
    end if
    
    c = b
    fc = fb
    do iIteration = 1, RF % MaximumIteration 
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
      call RF % FunctionEvaluator ( RF % FunctionParameters, b, fb )
    end do
    
    RF % SolutionAccuracy = abs ( fb )
    Root=b
    
    RF % Success = .false.

  end subroutine SolveBrent
  
  
  subroutine SolveSecant ( RF, Guess_1, Guess_2, Root )
  
    class ( RootFindingForm ), intent ( inout ) :: &
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
    
    call RF % FunctionEvaluator ( RF % FunctionParameters, X_0, Y_0 )
    call RF % FunctionEvaluator ( RF % FunctionParameters, X_1, Y_1 )
    
    do iIteration = 1, RF % MaximumIteration
      
      RF % nIterations = iIteration
      
      X = X_1 - Y_1 * ( X_1 - X_0 ) / ( Y_1 - Y_0 )
      
      X_0 = X_1
      Y_0 = Y_1
      
      X_1 = X
      call RF % FunctionEvaluator ( RF % FunctionParameters, X_1, Y_1 )
      
      RF % SolutionAccuracy = abs ( Y_1 )
      
      if ( abs ( ( X_1 - X_0 ) /  X_0 ) <= RF % RequestedAccuracy &
           .or. abs ( Y_1 ) <= RF % RequestedAccuracy ) &
      then
        RF % Success = .true. 
        Root = X_1
        exit
      end if
      
    end do
    
  end subroutine SolveSecant
  
  
  Subroutine Finalize ( RF )
  
    type ( RootFindingForm ), intent ( inout ) :: &
      RF
    
    nullify ( RF % FunctionParameters )
    nullify ( RF % FunctionDerivativeEvaluator )
    nullify ( RF % FunctionEvaluator )
  
  end subroutine Finalize
  

end module RootFinding_Form
