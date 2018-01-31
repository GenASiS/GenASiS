module Jacobi_Form
  
  use ISO_C_BINDING
  use Basics
  
  implicit none
  private
  
  type, public ::  JacobiForm
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Input, &
      Output
    type ( c_ptr ) :: &
      D_Input, &
      D_Output
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute_GPU
    procedure, public, pass :: &
      Compute_CPU
    final :: &
      Finalize
  end type JacobiForm
  

contains

  subroutine Initialize ( J )
  
    class ( JacobiForm ), intent ( inout ), target :: &
      J
    
    integer ( KDI ), dimension ( 2 ) :: &
      nCells
    real ( KDR ), dimension ( : , : ), pointer :: &
      A_O
    type ( c_ptr ) :: &
      D_A_O
    
    nCells = [ 512, 512 ]
    call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
    
    allocate ( J % Input ( nCells ( 1 ), nCells ( 2 ) ) )
    allocate ( J % Output ( nCells ( 1 ), nCells ( 2 ) ) )
    
    A_O => J % Input
  
    call random_number ( J % Input )
    
    J % D_Input = C_NULL_PTR
    
    
  end subroutine Initialize
  
  
  subroutine Compute_GPU ( J, A_O, A_U )
  
    class ( JacobiForm ), intent ( inout ) :: &
      J
    real ( KDR ), dimension ( :, : ), intent ( in ), target  :: &
      A_O
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      A_U
    
    integer ( KDI ) :: &
      iI, &
      iV, jV, &
      nIterations
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    type ( c_ptr ) :: &
      D_A_O
    
    nIterations = 1000
    nV = shape ( A_O )
    
    !$OMP target enter data map ( to: A_O )
    !$OMP target enter data map ( alloc: A_U )
    
    !-- !$OMP target data map ( alloc : A_U ) map ( A_O ) 
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
    
    !-- !$OMP end target data
    
    print*, 'A_O Loc Before', c_loc ( A_O )
    print*, 'D_Input Before', J % D_Input
    
    !call Show ( A_O ( 1, 1 ), 'A_O Before', CONSOLE % INFO_1 )
    
    !$OMP target data map ( A_O ) use_device_ptr ( A_O )
    print*, 'A_O Inside', c_loc ( A_O )
    J % D_Input = c_loc ( A_O )
    !$OMP end target data
    
    print*, 'A_O Loc After', c_loc ( A_O )
    print*, 'D_Input After', J % D_Input
    
    !call Show ( A_O ( 1, 1 ), 'A_O After', CONSOLE % INFO_1 )
    
  end subroutine Compute_GPU
  

  subroutine Compute_CPU ( J )
  
    class ( JacobiForm ), intent ( inout ), target :: &
      J
    
    integer ( KDI ) :: &
      iI, &
      iV, jV, &
      nIterations
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    real ( KDR ), dimension ( :, : ), pointer :: &
      A_O, &
      A_U
      
    A_O => J % Input
    A_U => J % Output
    
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
    
    nullify ( A_O, A_U )
    
  end subroutine Compute_CPU
  
  
  subroutine Finalize ( J ) 
  
    type ( JacobiForm ), intent ( inout ) :: &
      J
      
    deallocate ( J % Output ) 
    deallocate ( J % Input ) 
  
  end subroutine Finalize
  

end module Jacobi_Form


program Difference_Form_Test

  use Basics
  use Jacobi_Form

  implicit none
  
  integer ( KDI ), dimension ( 2 ) :: &
    nCells
  type ( TimerForm ) :: &
    T_GPU, &
    T_CPU
  type ( JacobiForm ) :: &
    J
  
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'JacobiTest' )
  
  call J % Initialize ( )
  
  call Show ( 'Calling Jacobi', CONSOLE % INFO_1 )
  
  call T_GPU % Initialize ( 'Jacobi GPU ', Level = 1 )
  call T_CPU % Initialize ( 'Jacobi CPU ', Level = 1 )
  
  !call T_CPU % Start ( )
  !call J % Compute_CPU ( )
  !call T_CPU % Stop ( )
  
  call T_GPU % Start ( )
  call J % Compute_GPU ( J % Input, J % Output )
  call T_GPU % Stop ( )
  
  call T_GPU % ShowTotal ( CONSOLE % INFO_1 )
  call T_CPU % ShowTotal ( CONSOLE % INFO_1 )
  
  call Show ( T_CPU % TotalTime / T_GPU % TotalTime, 'GPU SpeedUp Factor' )
  
  deallocate ( PROGRAM_HEADER )

end program Difference_Form_Test
