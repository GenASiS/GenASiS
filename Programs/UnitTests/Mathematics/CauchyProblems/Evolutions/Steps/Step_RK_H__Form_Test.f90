program Step_RK_H__Form_Test

  !-- Step_RK_Header__Form_Test

  use Basics
  use Steps

  implicit none

  real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
    A
  real ( KDR ), dimension ( 2 : 2 ) :: &
    C
  real ( KDR ), dimension ( 1 : 2 ) :: &
    B
  type ( Step_RK_H_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Step_RK_H__Form_Test', DimensionalityOption = '2D' )

  call Clear ( A )
  A ( 2, 1 ) = 1.0_KDR

  B ( 1 ) = 0.5_KDR
  B ( 2 ) = 0.5_KDR

  C ( 2 ) = 1.0_KDR
    
  allocate ( S )
  call S % Initialize ( A, B, C, nEquations = 1 )

  call S % Show ( )

  deallocate ( S )
  deallocate ( PROGRAM_HEADER )

end program Step_RK_H__Form_Test
