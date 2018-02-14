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
    type ( TimerForm ) :: &
      Timer_GPU, &
      Timer_CPU, &
      Timer_DataTransfer
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute_GPU
    procedure, public, pass :: &
      Compute_CPU
    procedure, public, pass :: &
      Validate
    final :: &
      Finalize
  end type JacobiForm
  

contains

  subroutine Initialize ( J )
  
    class ( JacobiForm ), intent ( inout ), target :: &
      J
    
    integer ( KDI ) :: &
      iV, jV, &
      Error
    integer ( KDI ), dimension ( 2 ) :: &
      nCells
    real ( KDR ), dimension ( :, : ), pointer :: &
      I, O
    
    nCells = [ 4, 4 ]
    call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
    
    call J % Timer_GPU % Initialize ( 'Jacobi GPU ', Level = 1 )
    call J % Timer_CPU % Initialize ( 'Jacobi CPU ', Level = 1 )
    call J % Timer_DataTransfer % Initialize ( 'Data Transfer ', Level = 1 )
    
    allocate ( J % Input ( nCells ( 1 ), nCells ( 2 ) ) )
    allocate ( J % Output ( nCells ( 1 ), nCells ( 2 ) ) )
    
    call random_number ( J % Input )
    J % Output = 0.0_KDR
    
    !-- Allocate memory on Device
    J % D_Input  = Allocate_D ( size ( J % Input ) )
    J % D_Output = Allocate_D ( size ( J % Output ) )
    
    call J % Timer_DataTransfer % Start ( )
    
    !-- Update device's J % Input 
    I => J % Input
    Error = AssociateTarget ( c_loc ( I ), J % D_Input, size ( I ), 0 )
    !$OMP target update to ( I )
    Error = DisassociateTarget ( c_loc ( I ) )
    
    !-- Update device's J % Output
    O => J % Output
    Error = AssociateTarget ( c_loc ( O ), J % D_Output, size ( O ), 0 )
    !$OMP target update to ( O )
    Error = DisassociateTarget ( c_loc ( O ) )
    
    call J % Timer_DataTransfer % Stop ( )
    
    nullify ( O )
    nullify ( I )
    
  end subroutine Initialize
  
  
  subroutine Compute_GPU ( J )
  
    class ( JacobiForm ), intent ( inout ), target :: &
      J
    
    integer ( KDI ) :: &
      iI, &
      iV, jV, &
      nIterations, &
      Error
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    real ( KDR ), dimension ( :, : ), pointer :: &
      A_U, &
      A_O
          
    A_O => J % Input
    A_U => J % Output
    
    nIterations = 1000
    nV = shape ( A_O )
    
    Error = AssociateTarget ( c_loc ( A_U ), J % D_Output, size ( A_U ), 0 )
    Error = AssociateTarget ( c_loc ( A_O ), J % D_Input, size ( A_O ), 0 )

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
    
    Error = DisassociateTarget ( c_loc ( A_U ) )
    Error = DisassociateTarget ( c_loc ( A_O ) )
    
    nullify ( A_U )
    nullify ( A_O )
    
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
    
  end subroutine Compute_CPU
  
  
  subroutine Validate ( J )
  
    class ( JacobiForm ), intent ( inout ), target :: &
      J
  
    !-- Validate computation done on GPU with results on CPU
    
    integer ( KDI ) :: &
      Error
    real ( KDR ) :: &
      L1_Error
    real ( KDR ), dimension ( :, : ), allocatable, target :: &
      O_Host
    
    allocate &
      ( O_Host ( size ( J % Input, dim = 1 ), size ( J % Input, dim = 2 ) ) )
    
    call J % Timer_DataTransfer % Start ( )
    
    !-- Get data on GPU to host variable O_Host so we can validate
    Error &
      = AssociateTarget &
          ( c_loc ( O_Host ), J % D_Output, size ( O_Host ), 0 )
    !$OMP target update from ( O_Host )
    Error = DisassociateTarget ( c_loc ( O_Host ) )
    
    call J % Timer_DataTransfer % Stop ( )
    
    L1_Error = sum ( abs ( O_Host - J % Output ) ) / sum ( J % Output ) 
    
    call Show ( L1_Error,   'L1_Error' )
    !call Show ( O_Host,     'O_Host' )
    !call Show ( J % Input,  'J % Input' )
    !call Show ( J % Output, 'J % Output' )
  
  end subroutine Validate
  
  
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
  type ( JacobiForm ) :: &
    J
  
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'JacobiTest' )
  
  call J % Initialize ( )
  
  associate &
    ( T_GPU => J % Timer_GPU, &
      T_CPU => J % Timer_CPU )
        
  call Show ( 'Calling Jacobi', CONSOLE % INFO_1 )
  
  
  call T_GPU % Start ( )
  call J % Compute_GPU ( )
  call T_GPU % Stop ( )
  
  call T_CPU % Start ( )
  call J % Compute_CPU ( )
  call T_CPU % Stop ( )
  
  call J % Validate ( ) 
  
  call T_GPU % ShowTotal ( CONSOLE % INFO_1 )
  call T_CPU % ShowTotal ( CONSOLE % INFO_1 )
  call J % Timer_DataTransfer % ShowTotal ( CONSOLE % INFO_1 )
  
  call Show ( T_CPU % TotalTime / T_GPU % TotalTime, 'GPU SpeedUp Factor' )
  
  end associate
  
  deallocate ( PROGRAM_HEADER )

end program Difference_Form_Test
