program LinearEquations_LAPACK_Form_Test
 
  use Basics

  integer ( KDI ) :: &
    nEquations, &
    nSolutions, &
    Info
  integer ( KDI ), dimension ( : ), allocatable :: &
    iaPivot
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Matrix
  real ( KDR ), dimension ( :, : ), allocatable, target :: &
    RightHandSide
  real ( KDR ), dimension ( :, : ), pointer :: &
    Solution

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'LinearEquations_LAPACK_Form_Test', &
           AppendDimensionalityOption = .false. )

  nEquations = 3
  nSolutions = 1

  allocate ( iaPivot ( nEquations ) )
  allocate ( RightHandSide ( nEquations, nSolutions ) )
  allocate ( Matrix ( nEquations, nEquations ) )
  Solution => RightHandSide

  !-- Wikipedia "System of Linear Equations"
  Matrix &
    = reshape ( [ 3., 2., -1., 2., -2., 0.5, -1., 4.0, -1. ], &
                [ nEquations, nEquations ] )
  RightHandSide &
    = reshape ( [ 1., -2., 0. ], [ nEquations, nSolutions ] )
  call Show ( Matrix, 'Matrix' )
  call Show ( RightHandSide, 'RightHandSide' )

  call DGESV ( nEquations, nSolutions, Matrix, nEquations, iaPivot, &
               RightHandSide, nEquations, Info )

  call Show ( Solution, 'Solution' )
  call Show ( [ 1.0_KDR, -2.0_KDR, -2.0_KDR ], 'Expected Solution' )

  deallocate ( PROGRAM_HEADER )

end program LinearEquations_LAPACK_Form_Test
