module SelectDevice_Command
  
  use OMP_LIB
  use iso_c_binding
  use Specifiers
  use Device_C
  use OffloadEnabled_Function
  
  implicit none
  private
  
  public :: &
    SelectDevice
  
contains


  subroutine SelectDevice ( iDevice, ErrorOption )
  
    integer ( KDI ), intent ( in ) :: &
      iDevice
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    integer ( KDI ) :: &
      Error
    
    Error = -1
    
    if ( OffloadEnabled ( ) ) then
      call OMP_SET_DEFAULT_DEVICE ( iDevice )
      Error = SetDevice ( iDevice )
    end if
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine SelectDevice
  

end module SelectDevice_Command
