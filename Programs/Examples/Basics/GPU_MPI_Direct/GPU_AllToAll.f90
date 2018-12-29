program GPU_AllToAll

  !-- Test direct GPU to GPU communications
  !   FIXME: Only works with CCE so far
  !   RUNNING: export MPICH_RDMA_ENABLED_CUDA=1, use 1 MPI per device (node)

  use MPI
  use Basics
  
  implicit none
  
  integer ( KDI ) :: &
    iV
  real ( KDR ), dimension ( : ), allocatable :: &
    SendBuffer, &
    RecvBuffer
  
  
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'GPU_AllToAll' )
  
  associate ( C => PROGRAM_HEADER % Communicator )
  
  allocate ( SendBuffer ( C % Size * 4 ) )
  allocate ( RecvBuffer ( C % Size * 4 ) )
  SendBuffer = [ ( ( C % Rank + 1 ) * 10.0_KDR + iV, &
                     iV = 1, size ( SendBuffer ) ) ]
  RecvBuffer = 0.0_KDR
  
  call Show ( SendBuffer, 'SendBuffer Init' )
  
  !$acc data copyin ( SendBuffer ), copyout ( RecvBuffer )
  
  !-- SendBuffer is copied to device, set the one on host to -huge()
  SendBuffer = - huge ( KDR ) !-- Set this on host
  
  !$acc host_data use_device ( SendBuffer, RecvBuffer ) 
  !-- 'SendBuffer' and 'RecvBuffer' now refer to the variable in device
  
  !-- The following MPI_ALLTOALL is direct GPU to GPU
  call MPI_ALLTOALL &
         ( SendBuffer, 4, MPI_DOUBLE_PRECISION, &
           RecvBuffer, 4, MPI_DOUBLE_PRECISION, &
           C % Handle, C % Error )
  
  !$acc end host_data
  
  !-- The following variables are on host and should display -huge() and zero
  call Show ( SendBuffer, 'SendBuffer CPU (-huge)' )
  call Show ( RecvBuffer, 'RecvBuffer CPU (zeroes)' )
  
  !$acc end data
  
  !-- The following variable has been populated by the value from GPU
  call Show ( RecvBuffer, 'RecvBuffer Final' )
  
  end associate
  
  deallocate ( PROGRAM_HEADER )

end program GPU_AllToAll
