module Jacobi_Command

  use Basics
  
  implicit none
  private
  
  public :: &
    Jacobi_GPU, &
    Jacobi_CPU

contains

  subroutine Jacobi_GPU ( A_O, A_U )
  
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      A_O
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      A_U

    integer ( KDI ) :: &
      iI, &
      iV, jV, &
      nIterations
    integer ( KDI ), dimension ( 2 ) :: &
      nV
      
    nIterations = 1000
    nV = shape ( A_O )
    
    !$OMP target data map ( alloc : A_U ) map ( A_O ) 
    do iI = 1, nIterations
      
      !$OMP target teams distribute parallel do collapse ( 2 ) schedule ( static, 1 )
      do jV = 2, nV ( 2 ) - 1 
        !-- !$OMP parallel do
        do iV = 2, nV ( 1 ) - 1
          A_U ( iV, jV ) & 
            = 0.25_KDR * (   A_O ( iV, jV - 1 ) + A_O ( iV, jV + 1 ) &
                           + A_O ( iV - 1, jV ) + A_O ( iV + 1, jV ) )
        end do
        !-- !$OMP end parallel do
      end do
      !$OMP end target teams distribute parallel do
    
    end do
    !$OMP end target data
    
  end subroutine Jacobi_GPU
  

  subroutine Jacobi_CPU ( A_O, A_U )
  
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      A_O
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      A_U

    integer ( KDI ) :: &
      iI, &
      iV, jV, &
      nIterations
    integer ( KDI ), dimension ( 2 ) :: &
      nV
      
    nIterations = 1000
    nV = shape ( A_O )
    
    do iI = 1, nIterations
      
      !$OMP parallel do
      do jV = 2, nV ( 2 ) - 1 
        do iV = 2, nV ( 1 ) - 1
          A_U ( iV, jV ) & 
            = 0.25_KDR * (   A_O ( iV, jV - 1 ) + A_O ( iV, jV + 1 ) &
                           + A_O ( iV - 1, jV ) + A_O ( iV + 1, jV ) )
        end do
      end do
      !$OMP end parallel do
    
    end do
    
  end subroutine Jacobi_CPU
  
  
end module Jacobi_Command


program Difference_Form_Test

  use Basics
  use Jacobi_Command

  implicit none
  
  integer ( KDI ), dimension ( 2 ) :: &
    nCells
  real ( KDR ), dimension ( :, : ), allocatable :: &
    A_O, &
    A_U
  type ( TimerForm ) :: &
    T_GPU, &
    T_CPU
  
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'JacobiTest' )
  
  nCells = [ 512, 512 ]
  call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
  
  allocate ( A_O ( nCells ( 1 ), nCells ( 2 ) ) )
  allocate ( A_U ( nCells ( 1 ), nCells ( 2 ) ) )
  
  call random_number ( A_O )
  
  call Show ( 'Calling Jacobi', CONSOLE % INFO_1 )
  
  call T_GPU % Initialize ( 'Jacobi GPU ', Level = 1 )
  call T_CPU % Initialize ( 'Jacobi CPU ', Level = 1 )
  
  
  call T_CPU % Start ( )
  call Jacobi_CPU ( A_O, A_U )
  call T_CPU % Stop ( )
  
  call T_GPU % Start ( )
  call Jacobi_GPU ( A_O, A_U )
  call T_GPU % Stop ( )
  
  call T_GPU % ShowTotal ( CONSOLE % INFO_1 )
  call T_CPU % ShowTotal ( CONSOLE % INFO_1 )
  
  call Show ( T_CPU % TotalTime / T_GPU % TotalTime, 'GPU SpeedUp Factor' )

  deallocate ( PROGRAM_HEADER )

end program Difference_Form_Test
