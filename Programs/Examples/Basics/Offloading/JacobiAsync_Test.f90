module Jacobi_Form
  
  use iso_c_binding
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
    type ( TimerForm ) :: &
      Timer_GPU, &
      Timer_CPU, &
      Timer_DataTransfer, &
      Timer_StartTransfer, &
      Timer_FinishTransfer
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
    
    nCells = [ 1024, 1024 ]
    call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
    
    call J % Timer_GPU % Initialize ( 'Jacobi GPU ', Level = 1 )
    call J % Timer_CPU % Initialize ( 'Jacobi CPU ', Level = 1 )
    !call J % Timer_DataTransfer % Initialize ( 'Data Transfer ', Level = 1 )
    call J % Timer_StartTransfer % Initialize ( 'Start Transfer ', Level = 1 )
    call J % Timer_FinishTransfer % Initialize ( 'Finish Transfer ', Level = 1 )
    
    allocate ( J % Input ( nCells ( 1 ), nCells ( 2 ) ) )
    allocate ( J % Output ( nCells ( 1 ), nCells ( 2 ) ) )
    
    call random_number ( J % Input )
    J % Output = 0.0_KDR
    
    !-- Allocate memory on Device
    call AllocateDevice ( J % Input, J % D_Input )
    call AllocateDevice ( J % Output, J % D_Output )
    
    call J % Timer_StartTransfer % Start ( )
    !-- Update device's J % Input 
    I => J % Input
    call AssociateDevice ( J % D_Input, I, ErrorOption = Error )
    !$OMP target update to ( I ) depend ( out: I ) nowait
    
    !-- Update device's J % Output
    O => J % Output
    call AssociateDevice ( J % D_Output, O, ErrorOption = Error )
    !$OMP target update to ( O ) depend ( out: O ) nowait
    call J % Timer_StartTransfer % Stop ( )
    
    
    call J % Timer_FinishTransfer % Start ( )
    !$OMP task depend ( in: O )
    call Show ( 'Finishing transfer I' )
    call DisassociateDevice ( I, Error )
    !$OMP end task
    
    !$OMP task depend ( in: I )
    call Show ( 'Finishing transfer O' )
    call DisassociateDevice ( O, Error )
    !$OMP end task
    call J % Timer_FinishTransfer % Stop ( )
    
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
    
    call AssociateDevice ( J % D_Output, A_U, ErrorOption = Error )
    call AssociateDevice ( J % D_Input, A_O, ErrorOption = Error )
    
!    do iI = 1, nIterations
!      !$OMP target teams distribute parallel do collapse ( 2 ) schedule ( static, 1 )
!      do jV = 2, nV ( 2 ) - 1
!        !-- !$OMP parallel do
!        do iV = 2, nV ( 1 ) - 1
!          A_U ( iV, jV ) & 
!            = 0.25_KDR * (   A_O ( iV, jV - 1 ) + A_O ( iV, jV + 1 ) &
!                           + A_O ( iV - 1, jV ) + A_O ( iV + 1, jV ) )
!        end do
!        !-- !$OMP end parallel do
!      end do
!      !$OMP end target teams distribute parallel do
!    end do
    
    call DisassociateDevice ( A_U, Error )
    call DisassociateDevice ( A_O, Error )
    
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
    
    !call J % Timer_DataTransfer % Start ( )
    
    !-- Get data on GPU to host variable O_Host so we can validate
    call AssociateDevice ( J % D_Output, O_Host, ErrorOption = Error )
    !$OMP target update from ( O_Host )
    call DisassociateDevice ( O_Host, Error )
    
    !call J % Timer_DataTransfer % Stop ( )
    
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
  !call J % Timer_DataTransfer % ShowTotal ( CONSOLE % INFO_1 )
  call J % Timer_StartTransfer % ShowTotal ( CONSOLE % INFO_1 )
  call J % Timer_FinishTransfer % ShowTotal ( CONSOLE % INFO_1 )
  
  call Show ( T_CPU % TotalTime / T_GPU % TotalTime, 'GPU SpeedUp Factor' )
  
  end associate
  
  deallocate ( PROGRAM_HEADER )

end program Difference_Form_Test
