module DeallocateDevice_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    DeallocateDevice
  
contains

  
  subroutine DeallocateDevice ( Device )
  
    type ( c_ptr ), intent ( inout ) :: &
      Device
    
    call DeallocateTarget ( Device )
  
  end subroutine DeallocateDevice


end module DeallocateDevice_Command
