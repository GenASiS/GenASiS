program Integral_Form_Test

  use Basics
  use Calculus

  implicit none

  real ( KDR ) :: &
    Expected, &
    Integral, &
    Parameters
  type ( IntegralForm ) :: &
    I

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Integral_Form_Test', AppendDimensionalityOption = .false. )
  
  Parameters  =  0.0_KDR
  call I % Initialize ( Parameters )

  !-- Integral 1

  I % Integrand  =>  Integrand_1  
  
  Expected  =  0.842700792949715_KDR
  call I % Compute ( 0.0_KDR, 1.0_KDR, Integral )

  call I % Show ( '*** Integral 1' )
  call Show ( Expected, 'Expected' )
  call Show ( Integral, 'Computed' )
  call Show ( abs ( ( Integral - Expected ) / Expected ), 'RelativeError' )

  !-- Integrand 2

  I % Integrand  =>  Integrand_2  
  
  Expected  =  CONSTANT % PI
  call I % Compute ( 0.0_KDR, 1.0_KDR, Integral )

  call I % Show ( '*** Integral 2' )
  call Show ( Expected, 'Expected' )
  call Show ( Integral, 'Computed' )
  call Show ( abs ( ( Integral - Expected ) / Expected ), 'RelativeError' )
       
  deallocate ( PROGRAM_HEADER )

contains


  function Integrand_1 ( Parameters, X ) result ( F )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      F

    real ( KDR ) :: &
      SqrtPi

    SqrtPi  =  sqrt ( CONSTANT % PI )

    F  =  ( 2.0_KDR / SqrtPi )  *  exp ( - X * X )

  end function Integrand_1


  function Integrand_2 ( Parameters, X ) result ( F )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      F

    F  =  4.0_KDR  *  sqrt ( 1 - X * X )

  end function Integrand_2


end program Integral_Form_Test
