program PROGRAM_HEADER_Singleton_Test

  use OMP_LIB
  use Specifiers
  use Display
  use PROGRAM_HEADER_Singleton

  implicit none
  
  integer ( KDI ) :: &
    Thread
  character ( LDF ) :: &
    Name = 'PROGRAM_HEADER_Singleton_Test'

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize &
         ( Name, AppendDimensionalityOption = .false. )
  
  call Show ( PROGRAM_HEADER % Communicator % Rank, 'Hello world from rank' )

  !$OMP parallel private ( Thread )
  Thread = OMP_GET_THREAD_NUM ( )
  call Show ( Thread, 'Hello world from thread' )
  !$OMP end parallel

  deallocate ( PROGRAM_HEADER )

end program PROGRAM_HEADER_Singleton_Test
