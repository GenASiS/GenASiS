#include "Preprocessor"
module DeviceAddress_Function
  
  !-- Return the device address of previously OpenMP-associated Host address

  use iso_c_binding
  use Specifiers
  use OnDevice_Function

  implicit none
  private

  interface DeviceAddress
    module procedure DeviceAddress_1D
  end interface DeviceAddress

  public :: &
    DeviceAddress

contains


  function DeviceAddress_1D ( Value ) result ( DA )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value

    type ( c_ptr ) :: &
      DA
      
    if ( OnDevice ( Value ) ) then
#ifdef ENABLE_OMP_OFFLOAD
      !$OMP target data use_device_ptr ( Value )
#endif
      DA = c_loc ( Value )
#ifdef ENABLE_OMP_OFFLOAD
      !$OMP end target data
#endif
    else 
      DA = c_null_ptr
    end if
      
  end function DeviceAddress_1D


end module DeviceAddress_Function
