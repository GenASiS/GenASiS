module UpdateDevice_Command
  
  use iso_c_binding
  use Specifiers
  use AssociateHost_Command
  use DisassociateHost_Command
  
  implicit none
  private
  
  public :: &
    UpdateDevice
  
  interface UpdateDevice
    module procedure UpdateDevice_KDR_1D
    module procedure UpdateDevice_KDR_2D
  end interface UpdateDevice
  
contains

  
  subroutine UpdateDevice_KDR_1D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateHost ( Device, Value, ErrorOption )
    !$OMP target update to ( Value ) 
    call DisassociateHost ( Value, ErrorOption )
  
  end subroutine UpdateDevice_KDR_1D


  subroutine UpdateDevice_KDR_2D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateHost ( Device, Value, ErrorOption )
    !$OMP target update to ( Value ) 
    call DisassociateHost ( Value, ErrorOption )
  
  end subroutine UpdateDevice_KDR_2D


end module UpdateDevice_Command
