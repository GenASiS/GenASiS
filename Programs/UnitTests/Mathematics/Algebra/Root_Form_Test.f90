program Root_Form_Test

  use Basics
  use Algebra
  
  implicit none
  
  integer ( KDI ) :: &
    iValue    
  real ( KDR ) :: &
    Expected, &
    Root
  type ( Real_1D_Form ) :: &
    Parameters
  type ( RootForm ) :: &
    R

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'RootFinder_Form_Test', AppendDimensionalityOption = .false. )
  
  call Parameters % Initialize ( 10 )
  Parameters % Value  =  [ ( acos ( -1.0_KDR ) * iValue, iValue = 1, 10 ) ]
  
  call R % Initialize ( Parameters )
  R % Zero  =>  Sine  
  
  Expected  =  2.0_KDR  *  CONSTANT % PI

  !-- solve with brent method
  call R % Solve &
              ( [ 1.5_KDR * acos ( - 1.0_KDR ), &
                  2.6_KDR * acos ( - 1.0_KDR ) ], Root )
  
  call R % Show ( '*** Brent Method' )
  call Show ( Expected, 'Expected' )
  call Show ( Root, 'Computed' )
  call Show ( abs ( ( Root - Expected ) / Expected ), 'RelativeError' )
  
  Root = huge ( 0.0_KDR )
  !-- solve with secant method
  call R % Solve &
              ( 1.5_KDR * acos ( - 1.0_KDR ), &
                2.6_KDR * acos ( - 1.0_KDR ), Root )
  
  call R % Show ( '*** Secant Method' )
  call Show ( Expected, 'Expected' )
  call Show ( Root, 'Computed' )
  call Show ( abs ( ( Root - Expected ) / Expected ), 'RelativeError' )

  Root = huge ( 0.0_KDR )
  !-- solve with newton-raphson method
  R % ZeroDerivative => Cosine
  call R % Solve &
              ( [ 1.5_KDR * acos ( - 1.0_KDR ), &
                  2.6_KDR * acos ( - 1.0_KDR ) ], Root )
  
  call R % Show ( '*** NewtonRaphson Method' )
  call Show ( Expected, 'Expected' )
  call Show ( Root, 'Computed' )
  call Show ( abs ( ( Root - Expected ) / Expected ), 'RelativeError' )

  deallocate ( PROGRAM_HEADER )


contains


  function Sine ( Parameters, X ) result ( F )

    class ( * ), intent ( inout ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      F
      
    select type ( P => Parameters )
    type is ( Real_1D_Form )
      call Show ( P % Value, 'Parameters', CONSOLE % INFO_7 )  
    end select
    
    F = sin ( X )
    
    call Show ( [ X, F ], 'sin Input - Result', CONSOLE % INFO_7 )
  
  end function Sine

  
  function Cosine ( Parameters, X ) result ( F )

    class ( * ), intent ( inout ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      F
      
    select type ( P => Parameters )
    type is ( Real_1D_Form )
      call Show ( P % Value, 'Parameters', CONSOLE % INFO_7 )  
    end select
    
    F = cos ( X )
    
    call Show ( [ X, F ], 'sin Input - Result', CONSOLE % INFO_7 )
  
  end function Cosine

  
end program Root_Form_Test
