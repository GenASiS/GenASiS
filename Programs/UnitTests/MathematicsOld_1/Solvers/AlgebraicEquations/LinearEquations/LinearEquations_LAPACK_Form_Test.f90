program LinearEquations_LAPACK_Form_Test
 
  use Basics
  use LinearEquations_LAPACK__Form

  type ( LinearEquations_LAPACK_Form ) :: &
    LE

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'LinearEquations_LAPACK_Form_Test', &
           AppendDimensionalityOption = .false. )

  call LE % Initialize ( nEquations = 3, nSolutions = 1 )

  !-- Wikipedia "System of Linear Equations"
  LE % Matrix &
    = reshape ( [ 3., 2., -1., 2., -2., 0.5, -1., 4.0, -1. ], &
                [ LE % nEquations, LE % nEquations ] )
  LE % RightHandSide &
    = reshape ( [ 1., -2., 0. ], [ LE % nEquations, LE % nSolutions ] )
  call Show ( LE % Matrix, 'Matrix' )
  call Show ( LE % RightHandSide, 'RightHandSide' )

  call LE % Solve ( )

  call Show ( LE % Solution, 'Solution' )
  call Show ( [ 1.0_KDR, -2.0_KDR, -2.0_KDR ], 'Expected Solution' )

  deallocate ( PROGRAM_HEADER )

end program LinearEquations_LAPACK_Form_Test
