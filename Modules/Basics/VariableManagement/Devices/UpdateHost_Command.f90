module UpdateHost_Command
  
  use iso_c_binding
  use Specifiers
  use AssociateDevice_Command
  use DisassociateDevice_Command
  
  implicit none
  private
  
  public :: &
    UpdateHost
  
  interface UpdateHost
    module procedure UpdateHost_KDR_1D
    module procedure UpdateHost_KDR_2D
  end interface UpdateHost
  
contains

  
  subroutine UpdateHost_KDR_1D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateDevice ( Device, Value, ErrorOption )
    !$OMP target update from ( Value )
    call DisassociateDevice ( Value, ErrorOption )
  
  end subroutine UpdateHost_KDR_1D


  subroutine UpdateHost_KDR_2D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateDevice ( Device, Value, ErrorOption )
    !$OMP target update from ( Value )
    call DisassociateDevice ( Value, ErrorOption )
  
  end subroutine UpdateHost_KDR_2D


end module UpdateHost_Command
