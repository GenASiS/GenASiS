module AllocateDevice_Command_Test_Form

  use iso_c_binding
  use Specifiers
  
  implicit none
  private
  
  interface Show
    module procedure ShowAddress_1D
  end interface Show
  
  public :: &
    Show

contains


  subroutine ShowAddress_1D ( Address, Description )

    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      Address
    character ( * ), intent ( in ) :: &
      Description
    
    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel
    character ( LDB ) :: &
      Buffer

    print '(a35)', trim ( Description )
    do i = 1, size ( Address )
      
      write ( IndexLabel, fmt = '( i7 )' ) i
      write ( Buffer, fmt = ' ( z64 ) ' ) Address ( i )
      Buffer = '0x' //  adjustl ( Buffer )
      print &
        '(a38, a32)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) = ', Buffer
    end do

  end subroutine ShowAddress_1D


end module AllocateDevice_Command_Test_Form


program AllocateDevice_Command_Test

  use iso_c_binding
  use Specifiers
  use AllocateDevice_Command
  use AllocateDevice_Command_Test_Form
  
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

end program AllocateDevice_Command_Test
