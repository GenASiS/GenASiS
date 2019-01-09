module AllocateHost_Command_Test_Form

  use iso_c_binding
  use Specifiers
  
  implicit none
  private
  
  interface Show
    module procedure Show_C_Pointer_1D
  end interface Show
  
  public :: &
    Show

contains


  subroutine Show_C_Pointer_1D ( C_Pointer, Description )

    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      C_Pointer
    character ( * ), intent ( in ) :: &
      Description
    
    integer ( KDI ) :: &
      i
    integer ( KBI ) :: &
      Address
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel
    character ( LDB ) :: &
      Buffer

    print '(a35)', trim ( Description )
    do i = 1, size ( C_Pointer )
      
      write ( IndexLabel, fmt = '( i7 )' ) i
      Address = transfer ( C_Pointer ( i ), 1_KBI )
      write ( Buffer, fmt = ' ( z64 ) ' ) Address
      Buffer = '0x' //  adjustl ( Buffer )
      print &
        '(a38, a32)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) = ', Buffer
    end do

  end subroutine Show_C_Pointer_1D


end module AllocateHost_Command_Test_Form


program AllocateHost_Command_Test

  use OMP_LIB
  use iso_c_binding
  use Specifiers
  use AllocateDevice_Command
  use DeallocateDevice_Command
  use UpdateDevice_Command
  use UpdateHost_Command
  use AllocateHost_Command
  use DeallocateHost_Command
  use AllocateHost_Command_Test_Form
  
  implicit none
  
  integer ( KDI ) :: &
    iV, &
    iVrbl, &
    nValues = 2000, &
    nVrbls
  real ( KDR ) :: &
    StartTime, &
    TotalTime, &
    DataSize_GB
  real ( KDR ), dimension ( : ), pointer :: &
    Pi_Value_1D  !-- pinned
  real ( KDR ), dimension ( :, : ), pointer :: &
    Pi_Value_2D  !-- pinned
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Pa_Value_2D  !-- pageable
  type ( c_ptr ) :: &
    d_Pi_Value_2D, & 
    d_Pa_Value_2D
  type ( c_ptr ), dimension ( : ), allocatable :: &
    d_Pi_Value_1D, & 
    d_Pa_Value_1D
  
  call AllocateHost ( Pi_Value_1D, [ 10 ] )
  call AllocateHost ( Pi_Value_2D, [ 5, 10 ] )
  
  Pi_Value_1D = [ ( iV * 2, iV = 1, size ( Pi_Value_1D ) ) ]
  Pi_Value_2D = reshape ( [ ( iV * 3, iV = 1, size ( Pi_Value_2D ) ) ], [ 5, 10 ] )
  
  !print*, 'Value_1D', Pi_Value_1D
  !print*, 'Value_2D', Pi_Value_2D
  
  call DeallocateHost ( Pi_Value_2D )
  call DeallocateHost ( Pi_Value_1D )
  
  call AllocateHost ( Pi_Value_2D, [ nValues, nValues ] )
  allocate ( Pa_Value_2D ( nValues, nValues ) )
  
  call random_number ( Pi_Value_2D )
  call random_number ( Pa_Value_2D )
  
  !--  Contiguous memory
  call AllocateDevice ( Pi_Value_2D, d_Pi_Value_2D )
  call AllocateDevice ( Pa_Value_2D, d_Pa_Value_2D )
  
  !-- Multiple 1D memory on device
  nVrbls = size ( Pi_Value_2D, dim = 2 )
  allocate ( d_Pi_Value_1D ( nVrbls ) )
  do iVrbl = 1, nVrbls
    call AllocateDevice ( Pi_Value_2D ( :, iVrbl ), d_Pi_Value_1D ( iVrbl ) )
  end do
  
  DataSize_GB = 1.0_KDR * nValues ** 2 * 8 / 1.0e9_KDR
  print*, 'Data size (GB)', DataSize_GB
   
  !-- Timing pinned memory update Host-to-Device
  StartTime = OMP_GET_WTIME ( )
  call UpdateDevice ( Pi_Value_2D, d_Pi_Value_2D )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Pinned memory Data Transfer'
  print*, 'Host-to-Device Time (s)        :', TotalTime
  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  !-- Timing pageable memory update Host-to-Device
!  StartTime = OMP_GET_WTIME ( )
!  call UpdateDevice ( Pa_Value_2D, d_Pa_Value_2D )
!  TotalTime = OMP_GET_WTIME ( ) - StartTime
!  print*, 'Pageable memory Data Transfer'
!  print*, 'Host-to-Device Time (s)        :', TotalTime
!  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime
!  print*, ''
  
  !-- Timing pageable memory update Host-to-Device loop
  StartTime = OMP_GET_WTIME ( )
  !$OMP parallel do
  do iVrbl = 1, nVrbls
    call UpdateDevice ( Pi_Value_2D ( :, iVrbl ), d_Pi_Value_1D ( iVrbl ) )
  end do
  !$OMP end parallel do
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Pinnned memory Data Transfer Loop'
  print*, 'Host-to-Device Time (s)        :', TotalTime
  print*, 'Host-to-Device Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  !-- Timing pinned memory update Device-To-Host
  StartTime = OMP_GET_WTIME ( )
  call UpdateHost ( d_Pi_Value_2D, Pi_Value_2D )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Pinned memory Data Transfer'
  print*, 'Device-to-Host Time (s)        :', TotalTime
  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  !-- Timing pageable memory update Device-To-Host
!  StartTime = OMP_GET_WTIME ( )
!  call UpdateHost ( d_Pa_Value_2D, Pa_Value_2D )
!  TotalTime = OMP_GET_WTIME ( ) - StartTime
!  print*, 'Pageable memory Data Transfer'
!  print*, 'Device-to-Host Time (s)        :', TotalTime
!  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime
!  print*, ''
  
  !-- Timing pageable memory update Device-To-Host Loop
  StartTime = OMP_GET_WTIME ( )
  !$OMP parallel do
  do iVrbl = 1, nVrbls
    call UpdateHost ( d_Pi_Value_1D ( iVrbl ), Pi_Value_2D ( :, iVrbl ) )
  end do
  !$OMP end parallel do
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Pinned memory Data Transfer Loop'
  print*, 'Device-to-Host Time (s)        :', TotalTime
  print*, 'Device-to-Host Bandwith (GB/s) :', DataSize_GB / TotalTime
  print*, ''
  
  call DeallocateDevice ( d_Pa_Value_2D )
  call DeallocateDevice ( d_Pi_Value_2D )
  do iVrbl = 1, nVrbls 
    call DeallocateDevice ( d_Pi_Value_1D ( iVrbl ) )
  end do
  
  
  call DeallocateHost ( Pi_Value_2D )
  deallocate ( Pa_Value_2D )
  
end program AllocateHost_Command_Test
