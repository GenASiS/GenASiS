module DeallocateDevice_Command_Test_Form

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


end module DeallocateDevice_Command_Test_Form


program DeallocateDevice_Command_Test

  use iso_c_binding
  use Specifiers
  use AllocateDevice_Command
  use DeallocateDevice_Command
  use DeallocateDevice_Command_Test_Form
  
  implicit none
  
  real ( KDR ), dimension ( 10 ) :: &
    Value_1D
  real ( KDR ), dimension ( 10, 5 ) :: &
    Value_2D
  character ( LDB ) :: &
    Buffer1, &
    Buffer2
  
  type ( c_ptr ) :: &
    dValue_1D, &
    dValue_2D
  
  
  dValue_1D = C_NULL_PTR
  dValue_2D = C_NULL_PTR
  
  call Show ( [ dValue_1D, dValue_2D ], 'Addresses pre allocation' )
  
  call AllocateDevice ( Value_1D, dValue_1D )
  call AllocateDevice ( Value_2D, dValue_2D )
  
  call Show ( [ dValue_1D, dValue_2D ], 'Addresses post allocation' )
  
  call DeallocateDevice ( dValue_1D )
  call DeallocateDevice ( dValue_2D )
  
  call Show ( [ dValue_1D, dValue_2D ], 'Addresses post deallocation' )

end program DeallocateDevice_Command_Test
