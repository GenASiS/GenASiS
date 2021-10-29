module SineFunction_Module

  use Basics

  implicit none
  private

  public :: SineFunction
  public :: CosineFunction

contains


  subroutine SineFunction ( Parameters, Input, Result )

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
  
  end subroutine SineFunction

  
  subroutine CosineFunction ( Parameters, Input, Result )

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
  
  end subroutine CosineFunction

  
end module SineFunction_Module


program RootFinder_Form_Test

  use Basics
  use Algebra
  use SineFunction_Module
  
  implicit none
  
  integer ( KDI ) :: &
    iValue    
  real ( KDR ) :: &
    Root
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
  RF % ZeroFunction => SineFunction  
  
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
  RF % ZeroFunctionDerivative => CosineFunction
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
