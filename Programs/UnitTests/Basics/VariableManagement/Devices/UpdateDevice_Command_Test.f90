program UpdateDevice_Command_Test

  use iso_c_binding
  use omp_lib
  use Specifiers
  use AllocateDevice_Command
  use UpdateDevice_Command
  
  implicit none
  
  integer ( KDI ) :: &
    Error_1, &
    Error_2
  real ( KDR ) :: &
    StartTime, &
    TotalTime
  real ( KDR ), dimension ( : ), allocatable :: &
    Value_1D
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Value_2D
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
  
  allocate ( Value_1D ( 100 ) )
  allocate ( Value_2D ( 100, 100 ) )
  call AllocateDevice ( Value_1D, dValue_1D )
  call AllocateDevice ( Value_2D, dValue_2D )
  
  call random_number ( Value_1D )
  call random_number ( Value_2D )
  
  StartTime = OMP_GET_WTIME ( )
  call UpdateDevice ( Value_1D, dValue_1D, ErrorOption = Error_1 )
  call UpdateDevice ( Value_2D, dValue_2D, ErrorOption = Error_2 )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Data Transfer Time (Sync):', TotalTime
  
  StartTime = OMP_GET_WTIME ( )
  call UpdateDeviceAsync ( Value_1D, dValue_1D, ErrorOption = Error_1 )
  call UpdateDeviceAsync ( Value_2D, dValue_2D, ErrorOption = Error_2 )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Data Transfer Time (Async)', TotalTime
  
  StartTime = OMP_GET_WTIME ( )
  call FinishUpdate ( dValue_1D)
  call FinishUpdate ( dValue_2D)
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Finish Transfer Time', TotalTime
  
  
end program UpdateDevice_Command_Test
