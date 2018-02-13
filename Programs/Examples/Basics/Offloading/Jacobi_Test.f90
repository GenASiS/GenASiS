module Jacobi_Form
  
  use ISO_C_BINDING
  use Basics
  use OMP_LIB
  
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
    
    integer ( KDI ) :: &
      iV, jV, Error
    real ( KDR ), dimension ( :, : ), pointer :: &
      A_O, A_U
    integer ( KDI ), dimension ( 2 ) :: &
      nCells
    real ( KDR ), dimension ( :, : ), allocatable, target :: &
      LocalArray
    
    nCells = [ 4, 4 ]
    call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
    
    allocate ( J % Input ( nCells ( 1 ), nCells ( 2 ) ) )
    allocate ( J % Output ( nCells ( 1 ), nCells ( 2 ) ) )
    
    J % Input  = +1.0_KDR
    J % Output = -2.0_KDR
    
    J % D_Input = C_NULL_PTR
    
    !associate ( A_O => J % Input, A_U => J % Output )
    A_O => J % Input
    A_U => J % Output
    
    !-- !$OMP target enter data map ( to: A_O )
    !-- !$OMP target enter data map ( alloc: A_U )
    
    print*, 'A_O Loc Before', c_loc ( A_O )
    print*, 'A_U Loc Before', c_loc ( A_U )
    print*, 'D_Input Before', J % D_Input
    print*, 'D_Output Before', J % D_Output
    
    call Show ( A_O ( 1, 1 ), 'A_O Before', CONSOLE % INFO_1 )
    
    !-- !$OMP target data map ( A_O ) use_device_ptr ( A_O )
    print*, 'A_O Inside', c_loc ( A_O )
    J % D_Input  = Allocate_D ( size ( LocalArray ) )
    !-- !$OMP end target data
    
    !-- !$OMP target data map ( A_U ) use_device_ptr ( A_U )
    !print*, 'A_U Inside', c_loc ( A_U )
    J % D_Output = Allocate_D ( size ( A_U ) )
    !-- !$OMP end target data
    
    print*, 'A_O Loc After', c_loc ( A_O )
    print*, 'D_Input After', J % D_Input
    print*, 'A_U Loc After', c_loc ( A_U )
    write ( *, fmt = '( a10, z20 )' ) 'D_Output After', J % D_Output
    
    call Show ( A_O ( 1, 1 ), 'A_O After', CONSOLE % INFO_1 )
    
    !end associate
    
    allocate ( LocalArray ( nCells ( 1 ), nCells ( 2 ) ) )
    LocalArray = 0.0_KDR
    Error = AssociateTarget &
               ( c_loc ( LocalArray ), J % D_Input, size ( LocalArray ), 0 )
    call Show ( Error, 'Error 1' )
    
    !$OMP target teams distribute parallel do
    do jV = 1, nCells ( 2 )
      do iV = 1, nCells ( 1 ) 
        LocalArray ( iV, jV ) = 1.11111_KDR
      end do
    end do
    !$OMP end target teams distribute parallel do
    
    Error = DisassociateTarget ( c_loc ( LocalArray ) )
    call Show ( LocalArray, 'LocalArray' )
    
  end subroutine Initialize
  
  
  subroutine Compute_GPU ( J, A_Old, A_Update )
  
    class ( JacobiForm ), intent ( inout ), target :: &
      J
    real ( KDR ), dimension ( :, : ), intent ( in ), target  :: &
      A_Old
    real ( KDR ), dimension ( :, : ), intent ( out ), target :: &
      A_Update
    
    integer ( KDI ) :: &
      iI, &
      iV, jV, &
      nIterations, &
      Error
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    type ( c_ptr ) :: &
      D_A_O
    real ( KDR ), dimension ( :, : ), pointer :: &
      A_O_P
      
    call Show ( A_Update, 'A_Update Begin' )
    call Show ( A_Old, 'A_Old Begin' )
    
    nIterations = 1000
    nV = shape ( A_Old )
    
    Error = AssociateTarget &
               ( c_loc ( A_Update ), J % D_Output, size ( A_Update ), 0 )
    call Show ( Error, 'Error 2' )
    
    Error = AssociateTarget &
               ( c_loc ( A_Old ), J % D_Input, size ( A_Old ), 0 )
    call Show ( Error, 'Error 3' )
    
    !call random_number ( A_Old ) !-- Should have no effect, 
                               !   since A_Old is already on device
    
    do iI = 1, nIterations
      !$OMP target teams distribute parallel do collapse ( 2 ) schedule ( static, 1 )
      do jV = 1, nV ( 2 ) - 1
        !-- !$OMP parallel do
        do iV = 2, nV ( 1 ) - 1
          A_Update ( iV, jV ) & 
            = 0.25_KDR * (   A_Old ( iV, jV - 1 ) + A_Old ( iV, jV + 1 ) &
                           + A_Old ( iV - 1, jV ) + A_Old ( iV + 1, jV ) )
        end do
        !-- !$OMP end parallel do
      end do
      !$OMP end target teams distribute parallel do
    end do
    
    call Show ( A_Update, 'A_Update After Compute' )
    !$OMP target update from ( A_Update )
    call Show ( A_Update, 'A_Update After Update' )
    
    Error = DisassociateTarget ( c_loc ( A_Update ) )
    Error = DisassociateTarget ( c_loc ( A_Old ) )
    
    
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
    
    !call Show ( A_U, 'A_U CPU' )
    
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
  
  call T_CPU % Start ( )
  call J % Compute_CPU ( )
  call T_CPU % Stop ( )
  
  J % Output = -2.0_KDR
  
  call T_GPU % Start ( )
  call J % Compute_GPU ( J % Input, J % Output )
  call T_GPU % Stop ( )
  
  call T_GPU % ShowTotal ( CONSOLE % INFO_1 )
  call T_CPU % ShowTotal ( CONSOLE % INFO_1 )
  
  call Show ( T_CPU % TotalTime / T_GPU % TotalTime, 'GPU SpeedUp Factor' )
  
  deallocate ( PROGRAM_HEADER )

end program Difference_Form_Test
