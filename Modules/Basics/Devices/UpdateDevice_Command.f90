module UpdateDevice_Command
  
  use iso_c_binding
  use omp_lib
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    UpdateDevice
  
  interface UpdateDevice
    module procedure UpdateDevice_KDI
    module procedure UpdateDevice_KDR_1D
    module procedure UpdateDevice_KDR_2D
    module procedure UpdateDevice_KDR_3D
    module procedure UpdateDevice_KDR_4D
    module procedure UpdateDevice_KDL
  end interface UpdateDevice
  
contains

  
  subroutine UpdateDevice_KDI ( Value, Device, ErrorOption )
  
    integer ( KDI ), dimension ( .. ), intent ( in ), target :: &
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
      
    Error = OMP_TARGET_MEMCPY &
              ( Device, c_loc ( Value ), c_sizeof ( Value ), &
                0_c_size_t, 0_c_size_t, OMP_GET_DEFAULT_DEVICE ( ), &
                OMP_GET_INITIAL_DEVICE ( ) )
    
  end subroutine UpdateDevice_KDI 
  
  
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


  subroutine UpdateDevice_KDR_3D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
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
  
  end subroutine UpdateDevice_KDR_3D
  
  
  subroutine UpdateDevice_KDR_4D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ), target :: &
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
  
  end subroutine UpdateDevice_KDR_4D
  
  
  subroutine UpdateDevice_KDL ( Value, Device, ErrorOption )
  
    logical ( KDL ), dimension ( .. ), intent ( in ), target :: &
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
      
    Error = OMP_TARGET_MEMCPY &
              ( Device, c_loc ( Value ), c_sizeof ( Value ), &
                0_c_size_t, 0_c_size_t, OMP_GET_DEFAULT_DEVICE ( ), &
                OMP_GET_INITIAL_DEVICE ( ) )
    
  end subroutine UpdateDevice_KDL


end module UpdateDevice_Command
