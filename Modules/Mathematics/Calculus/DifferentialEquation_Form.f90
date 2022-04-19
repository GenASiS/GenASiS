module DifferentialEquation_Form

  !-- After routines in Secs. 15.1 and 15.2 of Numerical Recipes 
  !   in Fortran (1989) 

  use Basics
   
  implicit none
  private
  
  type, public :: DifferentialEquationForm
    integer ( KDI ), private :: &
      IGNORABILITY
    integer ( KDI ) :: &
      MaxSteps
    real ( KDR ) :: &
      RequestedAccuracy, &
      SolutionAccuracy
    real ( KDR ), dimension ( : ), allocatable :: &
      Solution
    character ( LDL ) :: &
      SolverType
    procedure ( CS ), public, pointer, nopass :: &
      ComputeSlope => null ( )
    class ( * ), private, pointer :: &
      Parameters => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Integrate
    final :: &
      Finalize
  end type DifferentialEquationForm

  abstract interface 
    subroutine CS ( Parameters, X, Y, dYdX )
      use Basics
      implicit none
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        Y
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        dYdX
    end subroutine CS
  end interface

    private :: &
      Step_RK4

      private :: &
        StepRaw_RK4

    real ( KDR ), private, parameter :: &
      P_Grow    =  - 0.2_KDR, &
      P_Shrink  =  - 0.25_KDR, &
      F_Cor     =   1.0_KDR / 15.0_KDR, &
      Safety    =   0.9_KDR, &
      ErrorCon  =  ( 4.0_KDR / Safety ) ** ( 1 / P_Grow )

contains


  subroutine Initialize &
               ( DE, Parameters, nEquations, AccuracyOption, SolverTypeOption, &
                 MaxStepsOption, VerbosityOption )

    class ( DifferentialEquationForm ), intent ( inout ) :: &
      DE
    class ( * ), intent ( in ), target :: &
      Parameters
    integer ( KDI ) :: &
      nEquations
    real ( KDR ), intent ( in ), optional :: &
      AccuracyOption
    character ( * ), intent ( in ), optional :: &
      SolverTypeOption
    integer ( KDI ), intent ( in ), optional :: &
      MaxStepsOption, &
      VerbosityOption

    DE % IGNORABILITY  =  CONSOLE % INFO_3 
    if ( present ( VerbosityOption ) ) &
      DE % IGNORABILITY  =  VerbosityOption 

    DE % MaxSteps  =  100000
    if ( present ( MaxStepsOption ) ) &
      DE % MaxSteps  =  MaxStepsOption

    DE % RequestedAccuracy  =  epsilon ( 1.0_KDR ) * 1.0e4_KDR
    if ( present ( AccuracyOption ) ) &
      DE % RequestedAccuracy  =  AccuracyOption

    DE % SolverType  =  'RK4'
    if ( present ( SolverTypeOption ) ) &
      DE % SolverType  =  SolverTypeOption

    allocate ( DE % Solution ( nEquations ) )
    DE % Solution  =  0.0_KDR

    DE % Parameters  =>  Parameters
  
  end subroutine Initialize


  subroutine Integrate ( DE, X_Start, X_Finish, H_Start )

    class ( DifferentialEquationForm ), intent ( inout ) :: &
      DE
    real ( KDR ), intent ( in ) :: &
      X_Start, &
      X_Finish, &
      H_Start

    logical ( KDL ) :: &
      Reject
    integer ( KDI ) :: &
      nS, &   !-- number of steps
      nOk, &  !-- number of successful steps taken
      nBad    !-- number of unsuccessful steps taken
    real ( KDR ) :: &
      X, &
      H, &
      H_Out, &
      H_New, &
      Epsilon, &
      Error
    real ( KDR ), dimension ( size ( DE % Solution ) ) :: &
      Y, &
      Y_Scale, &
      dYdX

    X     =  X_Start
    H     =  H_Start
    Y     =  DE % Solution

    nOk   =  0
    nBad  =  0

    H_Out  =  0.0_KDR
    H_New  =  0.0_KDR

    Epsilon  =  DE % RequestedAccuracy

    do nS  =  1, DE % MaxSteps

      call DE % ComputeSlope ( DE % Parameters, X, Y, dYdX )
      if ( X  +  H  >  X_Finish ) &
         H  =  X_Finish - X
      select case ( trim ( DE % SolverType ) )
      case ( 'RK4' )
        Y_Scale  =  abs ( Y + abs ( H * dYdX ) )  +  sqrt ( tiny ( 0.0_KDR ) )
        call Step_RK4 ( DE, Y, dYdX, X, H, Epsilon, Y_Scale, H_Out, H_New )
      case default
        call Show ( 'SolverType not recognized', CONSOLE % ERROR )
        call Show ( 'DifferentialEquation_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Integrate', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select

      if ( H_Out  ==  H ) then 
        nOk  =  nOk + 1
      else
        nBad  =  nBad + 1
      end if
      
      if ( X  >=  X_Finish ) then
        DE % Solution  =  Y
        call Show ( 'Solution obtained in', DE % IGNORABILITY )
        call Show ( nS, 'total steps', DE % IGNORABILITY )
        call Show ( nOk, 'good steps', DE % IGNORABILITY )
        call Show ( nBad, 'altered steps', DE % IGNORABILITY )
        return
      end if

      H  =  H_New

    end do !-- iS
    
    call Show ( 'MaxSteps exceeded' )
    call Show ( DE % MaxSteps, 'MaxSteps', CONSOLE % ERROR )
    call Show ( 'DifferentialEquation_Form', 'module', CONSOLE % ERROR )
    call Show ( 'Integrate', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine Integrate


  subroutine Finalize ( DE )
  
    type ( DifferentialEquationForm ), intent ( inout ) :: &
      DE
    
    if ( allocated ( DE % Solution ) ) &
      deallocate ( DE % Solution )

    nullify ( DE % Parameters )
    nullify ( DE % ComputeSlope )

  end subroutine Finalize  


  recursive subroutine Step_RK4 &
              ( DE, Y, dYdX, X, H_In, Epsilon_In, Y_Scale, H_Out, H_New )

    class ( DifferentialEquationForm ), intent ( inout ) :: &
      DE
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
      X_S, &  !-- X_Save
      H, &
      HH
    real ( KDR ), dimension ( size ( Y ) ) :: &
         Y_S, &  !-- Y_Save
         Y_T, &  !-- Y_Temp
      dYdX_S, &  !-- dYdX_Save
      Epsilon, &
      ErrorMax 

       X_S  =     X
       Y_S  =     Y
    dYdX_S  =  dYdX

    H   =  H_In
    HH  =  H  /  2.0_KDR

    call StepRaw_RK4 ( DE, Y_S, dYdX_S, HH, X_S, Y_T ) 

    X  =  X_S  +  HH

    call DE % ComputeSlope &
           ( DE % Parameters, X, Y_T, dYdX )

    call StepRaw_RK4 ( DE, Y_T, dYdX, HH, X, Y )

    X  =  X_S  +  H
    
    call StepRaw_RK4 ( DE, Y_S, dYdX_S, H, X_S, Y_T ) 

    ErrorMax = 0.0_KDR

    Y_T  =  Y  -  Y_T

    ErrorMax  =  max ( ErrorMax, abs ( Y_T  /  Y_Scale ) ) 
    Epsilon   =   Epsilon_In  *  Y_Scale

    ErrorMax  =  ErrorMax  /  Epsilon

    if ( maxval ( ErrorMax )  >  1.0_KDR ) then
      H  =  Safety  *  H  *  ( maxval ( ErrorMax ) ** P_Shrink )
      call Step_RK4 &
             ( DE, Y_S, dYdX_S, X_S, H, Epsilon_In, Y_Scale, H_Out, H_New )
    else
      H_Out  =  H
      if ( maxval ( ErrorMax )  >  ErrorCon ) then
        H_New = Safety  *  H  *  ( maxval ( ErrorMax ) ** P_Grow )
      else
        H_New = 4.0_KDR * H
      end if
      Y  =  Y  +  Y_T  *  F_Cor
    end if

  end subroutine Step_RK4
  

  subroutine StepRaw_RK4 ( DE, Y, dYdX, H, X, Y_Out ) 

    type ( DifferentialEquationForm ), intent ( inout ) :: &
      DE    
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y, &
      dYdX
    real ( KDR ), intent ( in ) :: &
      H, &
      X
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Y_Out
    
    real ( KDR ) :: &
      HH, &
      H6, &
      XH
    real ( KDR ), dimension ( size ( Y ) ) :: &
       YT, &
      dYT, &
      dYM

    HH  =  H  /  2.0_KDR
    H6  =  H  /  6.0_KDR
    XH  =  X  +  HH

    !-- First step
    YT  =  Y  +  HH  *  dYdX

    !-- Second step
    call DE % ComputeSlope ( DE % Parameters, XH, YT, dYT )
    YT  =  Y  +  HH  *  dYT

    !-- Third Step
    call DE % ComputeSlope ( DE % Parameters, XH, YT, dYM )
     YT  =   Y   +  H  *  dYM
    dYM  =  dYT  +  dYM

    !-- Fourth Step
    call DE % ComputeSlope ( DE % Parameters, X + H, YT, dYT )

    !-- Accumulate increments with proper weights
    Y_Out  =  Y  +  H6  *  ( dYdX  +  dYT  +  2.0_KDR  *  dYM )

  end subroutine StepRaw_RK4

  
end module DifferentialEquation_Form
