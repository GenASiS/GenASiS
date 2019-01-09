module OffloadEnabled_Function
  
  use iso_c_binding
  use Specifiers
  use Device_C

  implicit none
  private
  
  public :: &
    OffloadEnabled
  
contains


  function OffloadEnabled ( ) result ( OE )
    
    logical ( KDL ) :: &
      OE
    
    OE = IsOffloadEnabled ( )
  
  end function OffloadEnabled


end module OffloadEnabled_Function
