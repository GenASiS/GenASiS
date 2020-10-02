program DeviceAddress_Function_Test

  use iso_c_binding
  use Specifiers
  use AllocateDevice_Command
  use AssociateHost_Command
  use DeviceAddress_Function
  
  implicit none
  
  integer :: &
    Error_1
  real ( KDR ), dimension ( 10 ) :: &
    Value_1D
  type ( c_ptr ) :: &
    dValue_1D, &
    dCopy_1D
    
  call AllocateDevice ( Value_1D, dValue_1D )
  call AssociateHost ( dValue_1D, Value_1D, ErrorOption = Error_1 )
  
  dCopy_1D = DeviceAddress ( Value_1D )
  
  if ( c_associated ( dValue_1D, dCopy_1D ) ) then
    print*, 'SUCCESS: DeviceAddress return correct address'
    call ShowPointer ( [ dValue_1D, dCopy_1D ], 'Original and DeviceAddress' )
  else
    print*, 'FAIL: DeviceAddress return wrong address' 
    call ShowPointer ( [ dValue_1D, dCopy_1D ], 'Original and DeviceAddress' )
  end if
  

contains

  subroutine ShowPointer ( C_Pointer, Description )

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

  end subroutine ShowPointer

  
end program
