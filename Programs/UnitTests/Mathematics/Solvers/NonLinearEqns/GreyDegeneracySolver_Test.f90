module MyFunction_Module

  use Basics

  implicit none
  private

  public :: FunctionWrapper
  public :: DerivativeFunctionWrapper

contains

  subroutine FunctionWrapper ( Parameters, Input, Result )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      Input
    real ( KDR ), intent ( out ) :: &
      Result

    real ( KDR ) :: &
      Fermi_2, Fermi_3, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta

    call dfermi &
             ( 2.0_KDR, Input, 0.0_KDR, Fermi_2, &
               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call dfermi &
             ( 3.0_KDR, Input, 0.0_KDR, Fermi_3, &
               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta)

    select type ( P => Parameters )
    type is ( real ( KDR  ) )
    
    
    Result = ( 2.0_KDR * CONSTANT % PI ** 2.0_KDR ) ** ( 1.0_KDR / 3.0_KDR ) &
             * Fermi_3 / ( Fermi_2 ** ( 4.0_KDR / 3.0_KDR ) ) - P

    call Show ( [ Input, Result ], 'Input - Result', CONSOLE % INFO_7 )

    end select
  
  end subroutine FunctionWrapper

  
  subroutine DerivativeFunctionWrapper ( Parameters, Input, Result )

    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      Input
    real ( KDR ), intent ( out ) :: &
      Result
      
    real ( KDR ) :: &
      Fermi_2, Fermi_3, &
      fdeta_2, fdeta_3, &
      fdeta2, fdtheta, &
      fdtheta2, fdetadtheta
    
    call dfermi &
             ( 2.0_KDR, input, 0.0_KDR, Fermi_2, &
               fdeta_2, fdtheta, fdeta2, fdtheta2, fdetadtheta)
    call dfermi &
             ( 3.0_KDR, input, 0.0_KDR, Fermi_3, &
               fdeta_3, fdtheta, fdeta2, fdtheta2, fdetadtheta)
    
    Result = ( 2.0_KDR * CONSTANT % PI ** 2.0_KDR ) ** ( 1.0_KDR / 3.0_KDR ) &
             * Fermi_2 ** ( - 4.0_KDR / 3.0_KDR ) &
             * ( fdeta_3 - 4.0_KDR / 3.0_KDR * Fermi_3 / Fermi_2 * fdeta_2)
  
  end subroutine DerivativeFunctionWrapper

end module MyFunction_Module

program GreyDegeneracySolver_Test

  use Basics
  use NonLinearEqns
  use MyFunction_Module

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
    Root, &
    T, &
    Eta, &
    J, &
    N
  procedure ( FunctionEvaluatorInterface ), pointer :: &
    FunctionEvaluator => null ( ), &
    FunctionDerivativeEvaluator => null ( )
  real ( KDR ) :: &
    Parameters, &
    Fermi_2, Fermi_3, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta
  type ( RootFindingForm ) :: &
    RF

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'GreyDegeneracySolver_Test', AppendDimensionalityOption = .false. )

  FunctionEvaluator => FunctionWrapper
  FunctionDerivativeEvaluator => DerivativeFunctionWrapper

  T = 0.1_KDR
  Eta = -1.0_KDR

  N =  CONSTANT % PI  ** ( - 2.0_KDR ) * exp ( Eta ) * T ** 3.0_KDR

  J =  3.0_KDR * CONSTANT % PI ** ( 2.0_KDR / 3.0_KDR ) &
       * exp ( - Eta / 3.0_KDR ) * N ** ( 4.0_KDR / 3.0_KDR )

  call Show ( N, 'N' )
  call Show ( J, 'J' )

  Parameters = J / N ** ( 4.0_KDR / 3.0_KDR )

  call RF % Initialize ( Parameters, FunctionEvaluator, 1.0e-10_KDR )

  !-- solve with brent method
  call RF % Solve ( [ 9.0_KDR, 11.0_KDR ], Root )
  
  if ( RF % Success ) then
   call dfermi &
          ( 2.0_KDR, Root, 0.0_KDR, Fermi_2, &
            fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call dfermi &
           ( 3.0_KDR, Root, 0.0_KDR, Fermi_3, &
             fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta)
    call Show ( 'Brent Method' )
    call Show ( Root, 'Eta' )
    call Show ( J / N * Fermi_2 / Fermi_3, 'T')
    call Show ( RF % nIterations, 'nIterations' )
  end if
  
  Root = huge ( 0.0_KDR )
  !-- solve with secant method
  call RF % Solve ( 9.0_KDR, 11.0_KDR, Root )
  
  if ( RF % Success ) then
    call dfermi &
          ( 2.0_KDR, Root, 0.0_KDR, Fermi_2, &
            fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call dfermi &
           ( 3.0_KDR, Root, 0.0_KDR, Fermi_3, &
             fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta) 
    call Show ( 'Secant Method' )
    call Show ( Root, 'Eta' )
    call Show ( J / N * Fermi_2 / Fermi_3, 'T')
    call Show ( RF % nIterations, 'nIterations' )
  end if

  Root = huge ( 0.0_KDR )
  !-- solve with newton-raphson method
  call RF % Solve &
              ( FunctionDerivativeEvaluator, &
                 [ 9.0_KDR, 11.0_KDR ], Root )
  
  if ( RF % Success ) then
    call dfermi &
          ( 2.0_KDR, Root, 0.0_KDR, Fermi_2, &
            fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call dfermi &
           ( 3.0_KDR, Root, 0.0_KDR, Fermi_3, &
             fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta)
    call Show ( 'Newton Raphson Method' )
    call Show ( Root, 'Eta' )
    call Show ( J / N * Fermi_2 / Fermi_3, 'T')
    call Show ( RF % nIterations, 'nIterations' )
  end if

  deallocate ( PROGRAM_HEADER )


end program GreyDegeneracySolver_Test
