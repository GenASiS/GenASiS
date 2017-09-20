module LinearEquations_LAPACK__Form

  use Basics

  implicit none
  private

  type, public :: LinearEquations_LAPACK_Form
    integer ( KDI ) :: &
      nEquations, &
      nSolutions, &
      Info
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPivot
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Matrix, &
      RightHandSide
    real ( KDR ), dimension ( :, : ), pointer :: &
      Solution
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
  end type LinearEquations_LAPACK_Form

contains

  
  subroutine Initialize ( LE, nEquations, nSolutions )

    class ( LinearEquations_LAPACK_Form ), intent ( inout ), target :: &
      LE
    integer ( KDI ), intent ( in ) :: &
      nEquations, &
      nSolutions

    LE % nEquations  =  nEquations
    LE % nSolutions  =  nSolutions

    allocate ( LE % iaPivot ( LE % nEquations ) )
    allocate ( LE % RightHandSide ( LE % nEquations, LE % nSolutions ) )
    allocate ( LE % Matrix ( LE % nEquations, LE % nEquations ) )
    LE % Solution => LE % RightHandSide

  end subroutine Initialize


  subroutine Solve ( LE )

    class ( LinearEquations_LAPACK_Form ), intent ( inout ), target :: &
      LE

    call DGESV ( LE % nEquations, LE % nSolutions, LE % Matrix, &
                 LE % nEquations, LE % iaPivot, LE % RightHandSide, &
                 LE % nEquations, LE % Info )

  end subroutine Solve


  impure elemental subroutine Finalize ( LE )

    type ( LinearEquations_LAPACK_Form ), intent ( inout ) :: &
      LE

    nullify ( LE % Solution )
    
    if ( allocated ( LE % Matrix ) ) &
      deallocate ( LE % Matrix )
    if ( allocated ( LE % RightHandSide ) ) &
      deallocate ( LE % RightHandSide )
    if ( allocated ( LE % iaPivot ) ) &
      deallocate ( LE % iaPivot )

  end subroutine Finalize


end module LinearEquations_LAPACK__Form
