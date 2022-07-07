program DeviceInterface_Test
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  integer ( c_size_t ) :: &
    Free, &
    Total
  integer ( c_int ) :: &
    Error
  integer ( KDI ) :: &
    nValues
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Value
  
  Error = DeviceMemGetInfo ( Free, Total )
  print*, '====== DeviceMemGetInfo ======'
  print*, '  Free   : ', Free
  print*, '  Total  : ', Total
  print*, ''
  
  allocate ( Value ( 256 ** 3, 30 ) )
  print*, 'Allocating array with shape', shape ( Value )
  print*, ''
  call random_number ( Value )
  !$OMP target enter data map ( to : Value ) 
  !$OMP barrier 
  
  Error = DeviceMemGetInfo ( Free, Total )
  print*, '====== DeviceMemGetInfo ======'
  print*, '  Free   : ', Free
  print*, '  Total  : ', Total
  print*, ''

end program DeviceInterface_Test
