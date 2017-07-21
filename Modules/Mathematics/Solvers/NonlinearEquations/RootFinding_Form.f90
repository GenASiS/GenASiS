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
    procedure ( MD_FunctionEvaluatorInterface ), private, pointer, nopass :: &
      MD_FunctionEvaluator => null ( )
    class ( * ), private, pointer :: &
      FunctionParameters => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SolveBrent
    procedure, private, pass :: &
      SolveSecant
    procedure, private, pass :: &
      SolveNewtonRaphson
    procedure, private, pass :: &
      SolveBroyden
    generic :: &
      Solve => SolveBrent, SolveSecant, SolveNewtonRaphson, SolveBroyden
    final :: &
      Finalize
  end type RootFindingForm

  private :: &
    ComputeInverseJacobian
  
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

  abstract interface
    subroutine MD_FunctionEvaluatorInterface ( Parameters, Input, Result )
      use Basics
      implicit none
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        Input
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        Result
    end subroutine MD_FunctionEvaluatorInterface
  end interface
  
  public :: &
    ZeroFunctionEvaluatorInterface
  
contains

!-- Move Accuracy & Parameters to solve calls.
  subroutine Initialize &
               ( RF, FunctionParameters, FunctionEvaluatorOption, &
                 MD_FunctionEvaluatorOption, AccuracyOption, &
                 MaximumIterationOption, VerbosityOption  )
    class ( RootFindingForm ), intent ( inout ) :: &
      RF
    class ( * ), intent ( in ), target :: &
      FunctionParameters
    procedure ( ZeroFunctionEvaluatorInterface ), &
                intent ( in ), pointer, optional :: &
      FunctionEvaluatorOption
    procedure ( MD_FunctionEvaluatorInterface ), &
                 intent ( in ), pointer, optional :: &
      MD_FunctionEvaluatorOption
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
    
    if ( present ( FunctionEvaluatorOption ) ) then
      RF % FunctionEvaluator => FunctionEvaluatorOption
    else if ( present ( MD_FunctionEvaluatorOption ) ) then
      RF % MD_FunctionEvaluator => MD_FunctionEvaluatorOption
    else
      call Show ( ' No Function Evaluator Specified ', CONSOLE % ERROR )
     ! call PROGRAM_HEADER % abort
    end if 
  
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
  
  
  subroutine SolveNewtonRaphson ( RF, FunctionDerivativeEvaluator, Interval, Root )
    
    !-- Newton-Raphson Method, based on Numerical Recipes in
    !   Fortran (1992) routine "rtsafe" 
  
    class ( RootFindingForm ), intent ( inout ) :: &
      RF
    procedure ( ZeroFunctionEvaluatorInterface ), intent ( in ), pointer :: &
      FunctionDerivativeEvaluator
    real(KDR), dimension ( 2 ), intent ( in ) :: &
      Interval
    real(KDR), intent ( out ) :: &
      Root
      
    integer ( KDI ) :: &
      iIteration
    real ( KDR ) :: &
      x1, x2, xacc, &
      df, dx, dxold, f, fh, fl, temp, xh, xl

    RF % FunctionDerivativeEvaluator => FunctionDerivativeEvaluator
    
    xacc = RF % RequestedAccuracy ! epsilon ( 1.0_KDR )  * 1.0e1_KDR
    
    x1 = Interval(1)
    x2 = Interval(2)

    !call Function(Parameters, x1, fl)
    !call Derivative(Parameters, x1, df)
    !call Function(Parameters, x2, fh)
    !call Derivative(Parameters, x2, df)

    call RF % FunctionEvaluator ( RF % FunctionParameters, x1, fl )
    call Show ( fl, 'fl' )
    call RF % FunctionDerivativeEvaluator ( RF % FunctionParameters, x1, df ) 
    call RF % FunctionEvaluator ( RF % FunctionParameters, x2, fh )
    call Show ( fh, 'fh' )
    call RF % FunctionDerivativeEvaluator ( RF % FunctionParameters, x2, df ) 
    
    if ( ( fl > 0.0_KDR .and. fh > 0.0_KDR ) &
         .or. ( fl < 0.0_KDR .and. fh < 0.0_KDR ) ) then
      call Show &
             ( 'Invalid interval was given as arguments', RF % IGNORABILITY )
      call Show ( 'RootFindingForm', 'Class', RF % IGNORABILITY )
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

    call RF % FunctionEvaluator ( RF % FunctionParameters, Root, f )
    call RF % FunctionDerivativeEvaluator ( RF % FunctionParameters, Root, df ) 
    
    do iIteration = 1, RF % MaximumIteration
      RF % nIterations = iIteration

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
      call RF % FunctionEvaluator ( RF % FunctionParameters, Root, f )
      call RF % FunctionDerivativeEvaluator ( RF % FunctionParameters, Root, df ) 
      if ( f < 0.0_KDR ) then
        xl = Root
      else
        xh = Root
      endif
    end do
    
    call Show('FindRoot exceeded maximum iterations', RF % IGNORABILITY )
    call Show( &
           'FindRoot could not find root to the specified accuracy', &
           RF % IGNORABILITY )

     RF % Success = .false.

  end subroutine SolveNewtonRaphson


  subroutine SolveBroyden ( RF, Guess_1, Guess_2, Root )
  
     class ( RootFindingForm ), intent ( inout ) :: &
       RF
     real ( KDR ), dimension ( : ), intent ( in ) :: &
       Guess_1, &
       Guess_2
     real ( KDR ), dimension ( : ), intent ( out ) :: &
       Root
      
     
    integer ( KDI ) :: &
      iIteration, &
      i, &
      j
    real ( KDR ) :: &
      IJ_Y, &
      dX_L2, &
      X_0_L2, &
      X_Error, &
      Y_Error
    real ( KDR ), dimension ( : ), allocatable :: &
      X, &
      X_0, X_1, dX, &
      Y_0, Y_1, dY
    real ( KDR ), dimension ( :, : ), allocatable :: &
      IJ_0, &
      IJ_1
      
    RF % Success = .false.
    RF % nIterations = 0
    
    allocate &
      ( X ( size ( Guess_1 ) ), &
        X_0 ( size ( Guess_1 ) ), &
        X_1 ( size ( Guess_1 ) ), &
        dX  ( size ( Guess_1 ) ), &
        Y_0 ( size ( Guess_1 ) ), &
        Y_1 ( size ( Guess_1 ) ), &
        dY  ( size ( Guess_1 ) ), &
        IJ_0 ( size ( Guess_1 ), size ( Guess_1 ) ), &
        IJ_1 ( size ( Guess_1 ), size ( Guess_1 ) ) )

    IJ_0 = 0.0_KDR  
    IJ_1 = 0.0_KDR

    do i = 1, size ( Guess_1 )
      IJ_0 ( i, i ) = 1.0_KDR
    end do !-- i

    X_0 = Guess_1
    X_1 = Guess_2
    
    call RF % MD_FunctionEvaluator ( RF % FunctionParameters, X_0, Y_0 )
    call RF % MD_FunctionEvaluator ( RF % FunctionParameters, X_1, Y_1 )

    do iIteration = 1, RF % MaximumIteration
      
      RF % nIterations = iIteration

      dX = X_1 - X_0
      dY = Y_1 - Y_0

      call ComputeInverseJacobian ( IJ_0, dX, dY, IJ_1 )
      
      do i = 1, size ( Guess_1 )
        IJ_Y = 0.0_KDR
        do j = 1, size ( Guess_1 )
          IJ_Y = IJ_Y + IJ_1 ( i, j ) * Y_1 ( j )
        end do !-- j
        X ( i ) = X_1 ( i ) - IJ_Y
      end do !-- i
      
      X_0 = X_1
      Y_0 = Y_1
      IJ_0 = IJ_1
      
      X_1 = X

      call RF % MD_FunctionEvaluator ( RF % FunctionParameters, X_1, Y_1 )

      X_Error = 0.0_KDR
      Y_Error = 0.0_KDR

      do i = 1, size ( Guess_1 )
        dX_L2  = dX_L2 + ( X_1 ( i ) - X_0 ( i ) ) ** 2.0_KDR 
        X_0_L2 = X_0_L2 + X_0 ( i ) ** 2.0_KDR
        Y_Error = Y_Error + Y_1 ( i ) ** 2.0_KDR
      end do !-- i

      X_Error = sqrt ( dX_L2 ) / ( 1.0_KDR + sqrt ( X_0_L2 ) )
      Y_Error = sqrt ( Y_Error )

      RF % SolutionAccuracy = Y_Error

      if ( X_Error <= RF % RequestedAccuracy &
           .or. Y_Error <= RF % RequestedAccuracy ) then
        RF % Success = .true. 
        Root = X_1
        exit
      end if
      
    end do

    deallocate ( X, X_0, X_1, dX, Y_0, Y_1, dY, IJ_0, IJ_1 )
    
  end subroutine SolveBroyden

  
  Subroutine Finalize ( RF )
  
    type ( RootFindingForm ), intent ( inout ) :: &
      RF
    
    nullify ( RF % FunctionParameters )
    nullify ( RF % FunctionDerivativeEvaluator )
    nullify ( RF % FunctionEvaluator )
  
  end subroutine Finalize


  subroutine ComputeInverseJacobian ( IJ_0, dX, dY, IJ_1 )
    !-- Rudimentary Application of Sherman-Morrison formula 
   real ( KDR ), dimension ( :, : ), intent ( in ) :: &
     IJ_0
   real ( KDR ), dimension ( : ), intent ( in ) :: &
     dX, &
     dY
   real ( KDR ), dimension ( :, : ), intent ( out ) :: &
     IJ_1 
   
   integer ( KDI ) :: &
     i, &
     j
   real ( KDR ) :: &
     dX_IJ_dY
   real ( KDR ), dimension ( size ( dX ) ) :: &
     dX_IJ, &
     IJ_dY  
    
   IJ_dY = 0.0_KDR
   dX_IJ = 0.0_KDR
   dX_IJ_dY = 0.0_KDR

    do i = 1, size ( dX )
      do j = 1, size ( dX )
        IJ_dY ( i ) = IJ_dY ( i ) + IJ_0 ( i, j ) * dY ( j )
        dX_IJ ( i ) = dX_IJ ( i ) + dX ( j ) * IJ_0 ( j, i )
      end do !-- j
      dX_IJ_dY = dX_IJ_dY + dX_IJ ( i ) * dY ( i )
    end do !-- i

    dX_IJ_dY = max ( dX_IJ_dY, tiny ( 0.0_KDR ) )

    do i = 1, size ( dX )
      do j = 1, size ( dX )
        IJ_1 ( i, j ) = IJ_0 ( i, j ) &
                        + ( dX ( i ) - IJ_dY ( i ) ) / dX_IJ_dY * dX_IJ ( j )
      end do !-- j
    end do !-- i

  end subroutine ComputeInverseJacobian

  

end module RootFinding_Form
