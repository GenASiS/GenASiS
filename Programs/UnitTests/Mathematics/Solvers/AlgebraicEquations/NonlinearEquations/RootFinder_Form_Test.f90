module MySinusoidalFunction_Module

  use Basics

  implicit none
  private

  public :: MySinusoidalFunction
  public :: MySinusoidalDerivativeFunction

contains


  subroutine MySinusoidalFunction ( Parameters, Input, Result )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      Input
    real ( KDR ), intent ( out ) :: &
      Result
      
    select type ( P => Parameters )
    type is ( Real_1D_Form )
      call Show ( P % Value, 'Parameters', CONSOLE % INFO_7 )  
    end select
    
    Result = sin ( Input )
    
    call Show ( [ Input, Result ], 'sin Input - Result', CONSOLE % INFO_7 )
  
  end subroutine MySinusoidalFunction

  
  subroutine MySinusoidalDerivativeFunction ( Parameters, Input, Result )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      Input
    real ( KDR ), intent ( out ) :: &
      Result
      
    select type ( P => Parameters )
    type is ( Real_1D_Form )
      call Show ( P % Value, 'Parameters', CONSOLE % INFO_7 )  
    end select
    
    Result = cos ( Input )
    
    call Show ( [ Input, Result ], 'sin Input - Result', CONSOLE % INFO_7 )
  
  end subroutine MySinusoidalDerivativeFunction

end module MySinusoidalFunction_Module


program RootFinder_Form_Test

  use Basics
  use NonlinearEquations
  use MySinusoidalFunction_Module
  
  implicit none
  
  interface 
    subroutine FunctionEvaluatorInterface ( Parameters, Input, Result )
      use Basics
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        Input
      real ( KDR ), intent ( out ) :: &
        Result
    end subroutine FunctionEvaluatorInterface
  end interface 
  
  integer ( KDI ) :: &
    iValue    
  real ( KDR ) :: &
    Root
  procedure ( FunctionEvaluatorInterface ), pointer :: &
    FunctionEvaluator => null ( ), &
    FunctionDerivativeEvaluator => null ( )
  type ( Real_1D_Form ) :: &
    Parameters
  type ( RootFinderForm ) :: &
    RF

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RootFinder_Form_Test', AppendDimensionalityOption = .false. )
  
  call Parameters % Initialize ( 10 )
  Parameters % Value = [ ( acos ( -1.0_KDR ) * iValue, iValue = 1, 10 ) ]
  
  call RF % Initialize ( Parameters )
  RF % EvaluateZero => MySinusoidalFunction  
  
  !-- solve with brent method
  call RF % Solve &
              ( [ 1.5_KDR * acos ( - 1.0_KDR ), &
                  2.6_KDR * acos ( - 1.0_KDR ) ], Root )
  
  if ( RF % Success ) then
    call Show ( 'Brent Method' )
    call Show ( Root, 'Root' )
    call Show ( RF % nIterations, 'nIterations' )
    call Show ( RF % SolutionAccuracy, 'SolutionAccuracy' )
  end if
  
  Root = huge ( 0.0_KDR )
  !-- solve with secant method
  call RF % Solve &
              ( 1.5_KDR * acos ( - 1.0_KDR ), &
                2.6_KDR * acos ( - 1.0_KDR ), Root )
  
  if ( RF % Success ) then
    call Show ( 'Secant Method' )
    call Show ( Root, 'Root' )
    call Show ( RF % nIterations, 'nIterations' )
    call Show ( RF % SolutionAccuracy, 'SolutionAccuracy' )
  end if

  Root = huge ( 0.0_KDR )
  !-- solve with newton-raphson method
  RF % EvaluateDerivative => MySinusoidalDerivativeFunction
  call RF % Solve &
              ( [ 1.5_KDR * acos ( - 1.0_KDR ), &
                  2.6_KDR * acos ( - 1.0_KDR ) ], Root )
  
  if ( RF % Success ) then
    call Show ( 'Newton Raphson Method' )
    call Show ( Root, 'Root' )
    call Show ( RF % nIterations, 'nIterations' )
    call Show ( RF % SolutionAccuracy, 'SolutionAccuracy' )
  end if

  deallocate ( PROGRAM_HEADER )

end program RootFinder_Form_Test
