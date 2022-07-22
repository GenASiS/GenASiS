module GetDeviceMemoryInfo_Command

  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private 
  
  public :: &
    GetDeviceMemoryInfo
    
contains
  

  subroutine GetDeviceMemoryInfo ( Total, Used, Free, ErrorOption )
  
    type ( QuantityForm ), intent ( out ) :: &
      Total, &
      Used, &
      Free
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    integer ( c_size_t ) :: &
      FreeBytes, &
      TotalBytes
    integer ( KDI ) :: &
      Error
    
    call Total % Initialize ( 'MB', 0.0_KDR )
    call Used  % Initialize ( 'MB', 0.0_KDR )
    call Free  % Initialize ( 'MB', 0.0_KDR )
    
    Error = DeviceMemGetInfo ( FreeBytes, TotalBytes )
    
    if ( Error == 0 ) then
      Free  % Number = FreeBytes / ( 1.0e6_KDR )
      Total % Number = TotalBytes / ( 1.0e6_KDR )
      Used  % Number = Total - Free
    end if
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine GetDeviceMemoryInfo
  

end module GetDeviceMemoryInfo_Command
