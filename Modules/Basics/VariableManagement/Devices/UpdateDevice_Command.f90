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
  
    real ( KDR ), dimension ( : ), target, intent ( in ) :: &
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
                                
    !Address = transfer ( Device, 1_KBI )
    !write ( Buffer, fmt = ' ( z64 ) ' ) Address
    !Buffer = '0x' //  adjustl ( Buffer )
    !print*, '1D Pointer: ', trim ( Buffer )
    
    !call AssociateHost ( Device, Value, ErrorOption )
    !print*, 'Done associating'
    !-- !$OMP target update to ( Value )
    !print*, 'Done transferring' 
    !call DisassociateHost ( Value, ErrorOption )
    !print*, 'Done Disassociation' 
    
    Error = HostToDeviceCopyDouble &
              ( c_loc ( Value ), Device, size ( Value ), 0, 0 )
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine UpdateDevice_KDR_1D


  subroutine UpdateDevice_KDR_2D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), target, intent ( in ) :: &
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
                                
    !Address = transfer ( Device, 1_KBI )
    !write ( Buffer, fmt = ' ( z64 ) ' ) Address
    !Buffer = '0x' //  adjustl ( Buffer )
    !print*, '2D Pointer: ', trim ( Buffer )
             
    !call AssociateHost ( Device, Value, ErrorOption )
    !!$OMP target update to ( Value ) 
    !call DisassociateHost ( Value, ErrorOption )
    
    
    !print*, 'Start HostToDeviceCopy'
    Error = HostToDeviceCopyDouble &
              ( c_loc ( Value ), Device, size ( Value ), 0, 0 )
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
    !print*, 'Done HostToDeviceCopy: ', Error
  
  
  end subroutine UpdateDevice_KDR_2D


end module UpdateDevice_Command
