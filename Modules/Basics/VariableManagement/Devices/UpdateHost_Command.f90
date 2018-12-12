module UpdateHost_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  use AssociateHost_Command
  use DisassociateHost_Command
  
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
    real ( KDR ), dimension ( : ), intent ( out ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    integer ( KDI ) :: &
      Error
      
    !call AssociateHost ( Device, Value, ErrorOption )
    !!$OMP target update from ( Value )
    !call DisassociateHost ( Value, ErrorOption )
        
    Error = DeviceToHostCopyDouble &
              ( Device, c_loc ( Value ), size ( Value ), 0, 0 )

    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine UpdateHost_KDR_1D


  subroutine UpdateHost_KDR_2D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( :, : ), intent ( out ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( KDI ) :: &
      Error
  
    !call AssociateHost ( Device, Value, ErrorOption )
    !!$OMP target update from ( Value )
    !call DisassociateHost ( Value, ErrorOption )
    
    Error = DeviceToHostCopyDouble &
              ( Device, c_loc ( Value ), size ( Value ), 0, 0 )

    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine UpdateHost_KDR_2D


end module UpdateHost_Command
