module MySinusoidalFunction_Module

  use Basics

  implicit none
  private

  public :: MySinusoidalFunction

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

end module MySinusoidalFunction_Module


program RootFinding_Form_Test

  use Basics
  use NonLinearEqns
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
    FunctionEvaluator => null ( )
  type ( Real_1D_Form ) :: &
    Parameters
  type ( RootFindingForm ) :: &
    RF

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RootFinding_Form_Test', AppendDimensionalityOption = .false. )
  
  FunctionEvaluator => MySinusoidalFunction
  
  call Parameters % Initialize ( 10 )
  Parameters % Value = [ ( acos ( -1.0_KDR ) * iValue, iValue = 1, 10 ) ]
  
  call RF % Initialize ( Parameters, FunctionEvaluator, 0.5_KDR )
  
  !-- solve with brent method
  call RF % Solve &
              ( [ 1.5_KDR * acos ( - 1.0_KDR ), &
                  2.5_KDR * acos ( - 1.0_KDR ) ], Root )
  
  if ( RF % Success ) then
    call Show ( 'Brent Method' )
    call Show ( Root, 'Root' )
    call Show ( RF % nIterations, 'nIterations' )
  end if
  
  Root = huge ( 0.0_KDR )
  !-- solve with secant method
  call RF % Solve &
              ( 1.5_KDR * acos ( - 1.0_KDR ), &
                2.5_KDR * acos ( - 1.0_KDR ), Root )
  
  if ( RF % Success ) then
    call Show ( 'Secant Method' )
    call Show ( Root, 'Root' )
    call Show ( RF % nIterations, 'nIterations' )
  end if

  deallocate ( PROGRAM_HEADER )

end program RootFinding_Form_Test
!contains 

  ! subroutine MySinusoidalFunction ( Parameters, Input, Result )

  !   use Basics

  !   class ( * ), intent ( in ) :: &
  !     Parameters
  !   real ( KDR ), intent ( in ) :: &
  !     Input
  !   real ( KDR ), intent ( out ) :: &
  !     Result
      
  !   select type ( P => Parameters )
  !   type is ( ArrayReal_1D_Form )
  !     call Show ( P % Value, 'Parameters', CONSOLE % INFO_7 )  
  !   end select
    
  !   Result = sin ( Input )
    
  !   call Show ( [ Input, Result ], 'sin Input - Result', CONSOLE % INFO_7 )
  
  ! end subroutine MySinusoidalFunction
  

!end program RootFinding_Form_Test
