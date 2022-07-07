module NumberOfDevices_Function
  
  use OMP_LIB
  use Specifiers

  implicit none
  private
  
  public :: &
    NumberOfDevices
  
  interface NumberOfDevices
    module procedure NumberOfDevices_OMP
  end interface NumberOfDevices
  
contains


  function NumberOfDevices_OMP ( ) result ( ND )
  
    integer ( KDI ) :: &
      ND
    
#ifdef ENABLE_OMP_OFFLOAD
    ND = OMP_GET_NUM_DEVICES ( ) 
#else
    ND = 0
#endif
  
  end function NumberOfDevices_OMP


end module NumberOfDevices_Function
