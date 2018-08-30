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

  use iso_c_binding
  use omp_lib
  use Specifiers
  use AllocateDevice_Command
  use UpdateDevice_Command
  use AllocateHost_Command
  use DeallocateHost_Command
  use AllocateHost_Command_Test_Form
  
  implicit none
  
  integer ( KDI ) :: &
    iV, &
    nValues = 10000
  real ( KDR ) :: &
    StartTime, &
    TotalTime
  real ( KDR ), dimension ( : ), pointer :: &
    Pi_Value_1D  !-- pinned
  real ( KDR ), dimension ( :, : ), pointer :: &
    Pi_Value_2D  !-- pinned
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Pa_Value_2D  !-- pageable
  type ( c_ptr ) :: &
    d_Pi_Value_2D, & 
    d_Pa_Value_2D
  
  call AllocateHost ( Pi_Value_1D, [ 10 ] )
  call AllocateHost ( Pi_Value_2D, [ 5, 10 ] )
  
  Pi_Value_1D = [ ( iV * 2, iV = 1, size ( Pi_Value_1D ) ) ]
  Pi_Value_2D = reshape ( [ ( iV * 3, iV = 1, size ( Pi_Value_2D ) ) ], [ 5, 10 ] )
  
  print*, 'Value_1D', Pi_Value_1D
  print*, 'Value_2D', Pi_Value_2D
  
  call DeallocateHost ( Pi_Value_2D )
  call DeallocateHost ( Pi_Value_1D )
  
  call AllocateHost ( Pi_Value_2D, [ nValues, nValues ] )
  allocate ( Pa_Value_2D ( nValues, nValues ) )
  
  call random_number ( Pi_Value_2D )
  call random_number ( Pa_Value_2D )
  
  call AllocateDevice ( Pi_Value_2D, d_Pi_Value_2D )
  call AllocateDevice ( Pa_Value_2D, d_Pa_Value_2D )
  
  !-- Timing pinned memory update
  StartTime = OMP_GET_WTIME ( )
  !call UpdateDevice ( Pi_Value_2D, d_Pi_Value_2D )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Data Transfer Time (Pinned)   :', TotalTime
  
  !-- Timing pageable memory update
  StartTime = OMP_GET_WTIME ( )
  call UpdateDevice ( Pa_Value_2D, d_Pa_Value_2D )
  TotalTime = OMP_GET_WTIME ( ) - StartTime
  print*, 'Data Transfer Time (Pageable) :', TotalTime

end program AllocateHost_Command_Test
