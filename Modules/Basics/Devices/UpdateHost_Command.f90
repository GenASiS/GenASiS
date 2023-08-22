module UpdateHost_Command
  
  use iso_c_binding
  use omp_lib
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    UpdateHost
  
  interface UpdateHost
    module procedure UpdateHost_KDI
    module procedure UpdateHost_KDR_1D
    module procedure UpdateHost_KDR_2D
    module procedure UpdateHost_KDR_3D
  end interface UpdateHost
  
contains

  
  subroutine UpdateHost_KDI ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), dimension ( .. ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    integer ( KDI ) :: &
      Error
    integer ( KBI ) :: &
      Address
    character ( LDB ) :: &
      Buffer
      
    Error = OMP_TARGET_MEMCPY &
              ( c_loc ( Value ), Device, c_sizeof ( Value ), &
                0_c_size_t, 0_c_size_t, OMP_GET_INITIAL_DEVICE ( ), &
                OMP_GET_DEFAULT_DEVICE ( ) )
    
  end subroutine UpdateHost_KDI 

  
  subroutine UpdateHost_KDR_1D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( : ), intent ( out ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    integer ( KDI ) :: &
      Error
      
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
  
    Error = DeviceToHostCopyDouble &
              ( Device, c_loc ( Value ), size ( Value ), 0, 0 )

    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine UpdateHost_KDR_2D


  subroutine UpdateHost_KDR_3D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( :, :, : ), intent ( out ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( KDI ) :: &
      Error
  
    Error = DeviceToHostCopyDouble &
              ( Device, c_loc ( Value ), size ( Value ), 0, 0 )

    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine UpdateHost_KDR_3D


end module UpdateHost_Command
