program UpdateDevice_Command_Test

  use iso_c_binding
  use omp_lib
  use Specifiers
  use AllocateDevice_Command
  use UpdateDevice_Command
  use UpdateHost_Command
  
  implicit none
  
  integer ( KDI ) :: &
    Error_1, &
    Error_2, &
    nValues
  real ( KDR ) :: &
    StartTime, &
    TotalTime
  real ( KDR ), dimension ( : ), allocatable :: &
    Value_1D, &
    Value_1D_Save
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Value_2D, &
    Value_2D_Save
  character ( LDB ) :: &
    Buffer1, &
    Buffer2
  type ( c_ptr ) :: &
    dValue_1D, &
    dValue_2D
    
  Error_1 = - huge ( 1 )
  Error_2 = - huge ( 1 )
  
  dValue_1D = C_NULL_PTR
  dValue_2D = C_NULL_PTR
  
  nValues = 30000
  
  allocate ( Value_1D ( nValues ) )
  allocate ( Value_2D ( nValues, nValues ) )
  allocate ( Value_1D_Save ( nValues ) )
  allocate ( Value_2D_Save ( nValues, nValues ) )
  
  call AllocateDevice ( Value_1D, dValue_1D )
  call AllocateDevice ( Value_2D, dValue_2D )
  
  call random_number ( Value_1D )
  call random_number ( Value_2D )
  Value_1D_Save = Value_1D
  Value_2D_Save = Value_2D
  
  
  !-- Testing sync update
  StartTime = OMP_GET_WTIME ( )
  call UpdateDevice ( Value_1D, dValue_1D )
  call UpdateDevice ( Value_2D, dValue_2D )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Data Transfer Time (Sync):', TotalTime
  
  Value_1D = 0.0_KDR
  Value_2D = 0.0_KDR
  
  call UpdateHost ( dValue_1D, Value_1D )
  call UpdateHost ( dValue_2D, Value_2D )
  
  !-- check error
  print*, 'Error 1D', &
    sum ( abs ( Value_1D_Save - Value_1D ) ) &
      / sum ( abs ( Value_1D_Save ) )
  print*, 'Error 2D', & 
    sum ( abs ( Value_2D_Save - Value_2D ) ) &
      / sum ( abs ( Value_2D_Save ) )
      
  
  call random_number ( Value_1D )
  call random_number ( Value_2D )
  Value_1D_Save = Value_1D
  Value_2D_Save = Value_2D
  
  StartTime = OMP_GET_WTIME ( )
  !$OMP parallel
  !$OMP single nowait
  call UpdateDeviceAsync ( Value_1D, dValue_1D )
  call UpdateDeviceAsync ( Value_2D, dValue_2D )
  !$OMP taskwait
  !$OMP end parallel
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Data Transfer Time (Async)', TotalTime
  
  StartTime = OMP_GET_WTIME ( )
  call FinishUpdate ( dValue_1D )
  call FinishUpdate ( dValue_2D )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Finish Transfer Time', TotalTime
  

  Value_1D = 0.0_KDR
  Value_2D = 0.0_KDR
  
  call UpdateHost ( dValue_1D, Value_1D )
  call UpdateHost ( dValue_2D, Value_2D )
  
  !-- check error
  print*, 'Error 1D', &
    sum ( abs ( Value_1D_Save - Value_1D ) ) &
      / sum ( abs ( Value_1D_Save ) )
  print*, 'Error 2D', & 
    sum ( abs ( Value_2D_Save - Value_2D ) ) &
      / sum ( abs ( Value_2D_Save ) )
  
  
end program UpdateDevice_Command_Test
