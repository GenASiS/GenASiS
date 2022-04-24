module Root_Form

  use Basics
   
  implicit none
  private
  
  type, public :: RootForm 
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaxIterations, &
      nIterations
    real ( KDR ) :: &
      RequestedAccuracy, &
      Accuracy, &
      AbsolutePrecision, &
      RelativePrecision
    logical ( KDL ) :: &
      Success = .false.
    procedure ( Z ), public, pointer, nopass :: &
      Zero => null ( ), &
      ZeroDerivative => null ( )
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
    procedure, public, pass :: &
      Show => Show_R
    final :: &
      Finalize
  end type RootForm
  
  abstract interface
    function Z ( Parameters, X ) result ( F )
      use Basics
      implicit none
      class ( * ), intent ( inout ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ) :: &
        F
    end function Z
  end interface

  
contains


  subroutine Initialize &
               ( R, Parameters, PrecisionOption, MaxIterationsOption, &
                 VerbosityOption )

    class ( RootForm ), intent ( inout ) :: &
      R
    class ( * ), intent ( in ), target :: &
      Parameters
    real ( KDR ), intent ( in ), optional :: &
      PrecisionOption
    integer ( KDI ), intent ( in ), optional :: &
      MaxIterationsOption, &
      VerbosityOption
      
    R % IGNORABILITY  =  CONSOLE % INFO_3 
    if ( present ( VerbosityOption ) ) &
      R % IGNORABILITY  =  VerbosityOption 
    
    R % MaxIterations  =  20
    if ( present ( MaxIterationsOption ) ) &
      R % MaxIterations  =  MaxIterationsOption

    R % RequestedAccuracy  =  epsilon ( 1.0_KDR ) * 10.0_KDR 
    if ( present ( PrecisionOption ) ) &
      R % RequestedAccuracy  =  PrecisionOption
    
    R % Parameters  =>  Parameters
  
  end subroutine Initialize
  
  
  subroutine SolveSecant ( R, Guess_1, Guess_2, Root )
  
    class ( RootForm ), intent ( inout ) :: &
      R
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
      
    R % Success = .false.
    R % nIterations = 0
      
    X_0 = Guess_1
    X_1 = Guess_2
    
    Y_0 = R % Zero ( R % Parameters, X_0 )
    Y_1 = R % Zero ( R % Parameters, X_1 )
    
    do iIteration = 1, R % MaxIterations
      
      R % nIterations = iIteration
      
      R % AbsolutePrecision = abs ( Y_1 )
      Root = X_1

      R % Accuracy           =  abs ( Y_1 )
      R % AbsolutePrecision  =  abs ( X_1 - X_0 )
      R % RelativePrecision  =  abs ( ( X_1 - X_0 ) )  &
                                /  max ( abs ( X_1 ), tiny ( 0.0_KDR ) )

      if ( R % Accuracy  <=  R % RequestedAccuracy &
           .or. R % AbsolutePrecision  <=  R % RequestedAccuracy   &
           .or. R % RelativePrecision  <=  R % RequestedAccuracy ) &
      then
        R % Success = .true. 
        exit
      end if
      
      if ( Y_1 == Y_0 ) &
        exit

      X = X_1 - Y_1 * ( X_1 - X_0 ) / ( Y_1 - Y_0 )
      
      X_0 = X_1
      Y_0 = Y_1
      
      X_1 = X
      Y_1 = R % Zero ( R % Parameters, X_1 )
      
    end do
    
  end subroutine SolveSecant
  

  subroutine Solve_B_NR ( R, Interval, Root )

    class ( RootForm ), intent ( inout ) :: &
      R
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      Interval
    real ( KDR ), intent ( out ) :: &
      Root
    
    if ( associated ( R % ZeroDerivative ) ) then
      call R % SolveNewtonRaphson ( Interval, Root )
    else
      call R % SolveBrent ( Interval, Root )
    end if
  end subroutine Solve_B_NR


  subroutine SolveBrent ( R, Interval, Root )
  
    !-- Van Wijngaarden-Dekker-Brent Method, based on Numerical Recipes in
    !   Fortran (1992) routine "zbrent"    
    
    class ( RootForm ), intent ( inout ) :: &
      R
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      Interval
    real ( KDR ), intent ( out ) :: &
      Root
      
    integer ( KDI ) :: &
      iIteration
    real ( KDR ) :: &
      a,b,c,d,e,fa,fb,fc,p,q,t,s,Accuracy1,xm
    
    R % Success = .false.
    R % nIterations = 0
  
    !-- FIXME: haven't set this in this routine
    R % RelativePrecision  =  0.0_KDR

    a = Interval ( 1 )
    b = Interval ( 2 )
    
    fa = R % Zero ( R % Parameters, a )
    fb = R % Zero ( R % Parameters, b )
    
    if ( ( fa > 0.0_KDR .and. fb > 0.0_KDR ) &
         .or. ( fa < 0.0_KDR .and. fb < 0.0_KDR ) ) &
    then
      call Show &
             ( 'Invalid interval was given as arguments', R % IGNORABILITY )
      call Show ( 'RootForm', 'Class', R % IGNORABILITY )
      call Show ( 'Solve', 'Method', R % IGNORABILITY )
      call Show ( Interval, 'Interval', R % IGNORABILITY )
      R % Success = .false.
      return
    end if
    
    c = b
    fc = fb
    do iIteration = 1, R % MaxIterations 
      R % nIterations = iIteration
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
      Accuracy1=2.0_KDR*epsilon(1.0_KDR)*abs(b)+0.5_KDR*R % RequestedAccuracy
      xm=0.5_KDR*(c-b)
      if (abs(xm) <= Accuracy1 .or. fb == 0.0_KDR) then
        Root=b
        R % Success = .true.
        R % AbsolutePrecision = abs ( xm )
        R % Accuracy = abs ( fb )
        return
      end if
      if (abs(e) >= Accuracy1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
          p=2.0_KDR*xm*s
          q=1.0_KDR-s
        else
          q=fa/fc
          t=fb/fc
          p=s*(2.0_KDR*xm*q*(q-t)-(b-a)*(t-1.0_KDR))
          q=(q-1.0_KDR)*(t-1.0_KDR)*(s-1.0_KDR)
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
      fb = R % Zero ( R % Parameters, b )
    end do
    
    R % AbsolutePrecision = abs ( xm )
    R % Accuracy = abs ( fb )
    Root=b
    
  end subroutine SolveBrent
  
  
  subroutine SolveNewtonRaphson ( R, Interval, Root )
    
    !-- Newton-Raphson Method, based on Numerical Recipes in
    !   Fortran (1992) routine "rtsafe" 
  
    class ( RootForm ), intent ( inout ) :: &
      R
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      Interval
    real ( KDR ), intent ( out ) :: &
      Root
      
    integer ( KDI ) :: &
      iIteration
    real ( KDR ) :: &
      x1, x2, xacc, &
      df, dx, dxold, f, fh, fl, temp, xh, xl

    R % Success = .true.
    R % nIterations = 0
    
    !-- FIXME: haven't set this in this routine
    R % RelativePrecision  =  0.0_KDR

    xacc = R % RequestedAccuracy
    
    x1 = Interval(1)
    x2 = Interval(2)

    fl = R % Zero ( R % Parameters, x1 )
    df = R % ZeroDerivative ( R % Parameters, x1 ) 
    fh = R % Zero ( R % Parameters, x2 )
    df = R % ZeroDerivative ( R % Parameters, x2 ) 
    
    if ( ( fl > 0.0_KDR .and. fh > 0.0_KDR ) &
         .or. ( fl < 0.0_KDR .and. fh < 0.0_KDR ) ) then
      call Show &
             ( 'Invalid interval was given as arguments', R % IGNORABILITY )
      call Show ( 'RootForm', 'Class', R % IGNORABILITY )
      call Show ( 'Solve', 'Method', R % IGNORABILITY )
      call Show ( Interval, 'Interval', R % IGNORABILITY )
      R % Success = .false.
      return
    end if
   
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
    
    f = R % Zero ( R % Parameters, Root )
    df = R % ZeroDerivative ( R % Parameters, Root ) 
    
    do iIteration = 1, R % MaxIterations
      R % nIterations = iIteration
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
      f = R % Zero ( R % Parameters, Root )
      R % AbsolutePrecision  =  abs ( dx )
      R % Accuracy = abs ( f )
      if ( abs ( dx ) < xacc ) return
      df = R % ZeroDerivative ( R % Parameters, Root ) 
      if ( f < 0.0_KDR ) then
        xl = Root
      else
        xh = Root
      endif
    end do
    
    call Show ( 'Root exceeded maximum iterations', R % IGNORABILITY )
    call Show ( 'Root could not find root to the requested precision', &
                R % IGNORABILITY )

     R % Success = .false.

  end subroutine SolveNewtonRaphson

  
  subroutine Show_R ( R, Name )

    class ( RootForm ), intent ( in ) :: &
      R
    character ( * ), intent ( in ) :: &
      Name

    call Show ( Name )
    call Show ( R % Success, 'Success' )
    call Show ( R % nIterations, 'nIterations' )
    call Show ( R % RequestedAccuracy, 'RequestedAccuracy' )
    call Show ( R % Accuracy, 'Accuracy' )
    call Show ( R % AbsolutePrecision, 'AbsolutePrecision' )
    call Show ( R % RelativePrecision, 'RelativePrecision' )

  end subroutine Show_R


  subroutine Finalize ( R )
  
    type ( RootForm ), intent ( inout ) :: &
      R
    
    nullify ( R % Parameters )
    nullify ( R % ZeroDerivative )
    nullify ( R % Zero )
  
  end subroutine Finalize
  

end module Root_Form
