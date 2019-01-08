module GetNumberOfDevices_Function
  
  use OMP_LIB
  use Specifiers

  implicit none
  private
  
  public :: &
    GetNumberOfDevices
  
  interface GetNumberOfDevices
    module procedure GetNumberOfDevices_OMP
  end interface GetNumberOfDevices
  
contains


  function GetNumberOfDevices_OMP ( ) result ( ND )
  
    integer ( KDI ) :: &
      ND
    
#ifdef ENABLE_OMP_OFFLOAD
    ND = OMP_GET_NUM_DEVICES ( ) 
#else
    ND = 0
#endif
  
  end function GetNumberOfDevices_OMP


end module GetNumberOfDevices_Function
