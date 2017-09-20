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
  use NonlinearEquations
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
    DensityFactor, &
    Root, &
    T, &
    Eta, &
    J, J_B, &
    N, N_B
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

  DensityFactor  &
    =  4.0_KDR * CONSTANT % PI  /  ( 2.0_KDR * CONSTANT % PI ) ** 3

  T = 1.0_KDR
  Eta = -10.0_KDR

  call DFERMI ( 2.0_KDR, Eta, 0.0_KDR, Fermi_2, &
                fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
  call DFERMI ( 3.0_KDR, Eta, 0.0_KDR, Fermi_3, &
                fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta ) 

  N  =  DensityFactor  *  T ** 3  *  Fermi_2
  J  =  DensityFactor  *  T ** 4  *  Fermi_3

  N_B  =  DensityFactor  *  T ** 3  *  2 * exp ( Eta )
  J_B  =  DensityFactor  *  T ** 4  *  6 * exp ( Eta )

  call Show ( N, 'N' )
  call Show ( J, 'J' )
  call Show ( N_B, 'N_B' )
  call Show ( J_B, 'J_B' )

  Parameters = J / N ** ( 4.0_KDR / 3.0_KDR )

  call RF % Initialize ( Parameters, FunctionEvaluator )

  !-- solve with brent method
  call RF % Solve ( [ -11.0_KDR, -9.0_KDR ], Root )
  
  if ( RF % Success ) then
    call Show ( 'Brent Method' )
    call Show ( Root, 'Eta' )
    call DFERMI ( 2.0_KDR, Root, 0.0_KDR, Fermi_2, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call DFERMI ( 3.0_KDR, Root, 0.0_KDR, Fermi_3, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta ) 
    call Show ( J / N * Fermi_2 / Fermi_3, 'T')
    call Show ( RF % nIterations, 'nIterations' )
    call Show ( RF % SolutionAccuracy, 'SolutionAccuracy' )
    call Show ( -3 * log ( CONSTANT % Pi ** ( - 2. / 3. )  / 3. &
                           *  J  /  N ** ( 4. / 3. ) ), 'Eta_B' )
  else
    call Show ( 'Brent Method FAIL' )
  end if
  
  Root = huge ( 0.0_KDR )
  !-- solve with secant method
  call RF % Solve ( -11.0_KDR, -9.0_KDR, Root )
  
  if ( RF % Success ) then
    call Show ( 'Secant Method' )
    call Show ( Root, 'Eta' )
    call DFERMI ( 2.0_KDR, Root, 0.0_KDR, Fermi_2, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call DFERMI ( 3.0_KDR, Root, 0.0_KDR, Fermi_3, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta ) 
    call Show ( J / N * Fermi_2 / Fermi_3, 'T')
    call Show ( RF % nIterations, 'nIterations' )
    call Show ( RF % SolutionAccuracy, 'SolutionAccuracy' )
    call Show ( -3 * log ( CONSTANT % Pi ** ( - 2. / 3. )  / 3. &
                           *  J  /  N ** ( 4. / 3. ) ), 'Eta_B' )
  end if

  Root = huge ( 0.0_KDR )
  !-- solve with newton-raphson method
  call RF % Solve &
              ( FunctionDerivativeEvaluator, &
                 [ -11.0_KDR, -9.0_KDR ], Root )
  
  if ( RF % Success ) then
    call Show ( 'Newton Raphson Method' )
    call Show ( Root, 'Eta' )
    call DFERMI ( 2.0_KDR, Root, 0.0_KDR, Fermi_2, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    call DFERMI ( 3.0_KDR, Root, 0.0_KDR, Fermi_3, &
                  fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta ) 
    call Show ( J / N * Fermi_2 / Fermi_3, 'T')
    call Show ( RF % nIterations, 'nIterations' )
    call Show ( RF % SolutionAccuracy, 'SolutionAccuracy' )
    call Show ( -3 * log ( CONSTANT % Pi ** ( - 2. / 3. )  / 3. &
                           *  J  /  N ** ( 4. / 3. ) ), 'Eta_B' )
  else
    call Show ( 'Newton Raphson FAIL' )
  end if

  deallocate ( PROGRAM_HEADER )


end program GreyDegeneracySolver_Test
