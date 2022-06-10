program DifferentialEquation_Form_Test

  use Basics
  use Calculus

  implicit none

  real ( KDR ) :: &
    Parameters
  type ( DifferentialEquationForm ) :: &
    DE_RK4

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'DifferentialEquation_Form_Test', &
           AppendDimensionalityOption = .false. )
  
  Parameters  =  0.0_KDR

  !-- RK4

  call DE_RK4 % Initialize &
         ( Parameters, nEquations = 2, &
           SolverTypeOption = 'RK4', &
           VerbosityOption = CONSOLE % INFO_1 )

  DE_RK4 % ComputeSlope  =>  ComputeSlope

  call Show ( 'Testing RK4' )
  call TestIntegration ( DE_RK4 )

  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestIntegration ( DE )

    class ( DifferentialEquationForm ), intent ( inout ) :: &
      DE

    associate &
      ( X    =>  DE % Solution ( 1 ), &
        V    =>  DE % Solution ( 2 ), &
        T_1  =>  0.0_KDR, &
        T_2  =>  CONSTANT % PI, &
        H_1  =>  1.0e-3_KDR ) !CONSTANT % PI / 2.0_KDR )

    DE % Solution ( 1 )  =  2.0_KDR
    DE % Solution ( 2 )  =  0.0_KDR

    call DE % Integrate &
           ( X_Start  = T_1, &
             X_Finish = T_2, &
             H_Start  = H_1 )

    call Show (   2.0_KDR * cos ( T_2 ), 'Expected X(Pi)' )
    call Show ( - 2.0_KDR * sin ( T_2 ), 'Expected V(Pi)' )

    call Show ( X, 'Computed X(Pi)' )
    call Show ( V, 'Computed V(Pi)' )

    end associate !-- X, etc.

  end subroutine TestIntegration


  subroutine ComputeSlope ( Parameters, X, Y, dYdX )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      dYdX

    dYdX ( 1 )  =    Y ( 2 )
    dYdX ( 2 )  =  - Y ( 1 )
  
  end subroutine ComputeSlope


end program DifferentialEquation_Form_Test
