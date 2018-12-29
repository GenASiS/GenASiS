module UpdateDevice_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
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
  
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( KDI ) :: &
      Error
    integer ( KBI ) :: &
      Address
    character ( LDB ) :: &
      Buffer
                                
    Error = HostToDeviceCopyDouble &
              ( c_loc ( Value ), Device, size ( Value ), 0, 0 )
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
      
  end subroutine UpdateDevice_KDR_1D


  subroutine UpdateDevice_KDR_2D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( KDI ) :: &
      Error
    integer ( KBI ) :: &
      Address
    character ( LDB ) :: &
      Buffer
                                
    Error = HostToDeviceCopyDouble &
              ( c_loc ( Value ), Device, size ( Value ), 0, 0 )
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine UpdateDevice_KDR_2D


end module UpdateDevice_Command
