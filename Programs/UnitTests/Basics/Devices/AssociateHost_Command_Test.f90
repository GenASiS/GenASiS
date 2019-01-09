module AssociateHost_Command_Test_Form

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
      Address = transfer ( C_Pointer, 1_KBI )
      write ( Buffer, fmt = ' ( z64 ) ' ) Address
      Buffer = '0x' //  adjustl ( Buffer )
      print &
        '(a38, a32)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) = ', Buffer
    end do

  end subroutine Show_C_Pointer_1D


end module AssociateHost_Command_Test_Form


program AssociateHost_Command_Test

  use iso_c_binding
  use Specifiers
  use AllocateDevice_Command
  use AssociateHost_Command
  use AssociateHost_Command_Test_Form
  
  implicit none
  
  integer ( KDI ) :: &
    Error_1, &
    Error_2
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
  
  Error_1 = - huge ( 1 )
  Error_2 = - huge ( 1 )
  
  dValue_1D = C_NULL_PTR
  dValue_2D = C_NULL_PTR
  
  call Show ( [ dValue_1D, dValue_2D ], 'Addresses pre allocation' )
  
  call AllocateDevice ( Value_1D, dValue_1D )
  call AllocateDevice ( Value_2D, dValue_2D )
  
  call Show ( [ dValue_1D, dValue_2D ], 'Addresses post allocation' )
  
  call AssociateHost ( dValue_1D, Value_1D, ErrorOption = Error_1 )
  call AssociateHost ( dValue_2D, Value_2D, ErrorOption = Error_2 )
  
  call Show ( [ dValue_1D, dValue_2D ], 'Addresses post association' )
  print*, 'Associate Error 1: ', Error_1
  print*, 'Associate Error 2: ', Error_2

end program AssociateHost_Command_Test
