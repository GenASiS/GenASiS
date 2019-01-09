module DeallocateHost_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    DeallocateHost
  
  interface DeallocateHost
    module procedure DeallocateHost_KDR_1D
    module procedure DeallocateHost_KDR_2D
  end interface DeallocateHost
  
contains

  
  subroutine DeallocateHost_KDR_1D ( Value )
  
    real ( KDR ), dimension ( : ), target, intent ( in ) :: &
      Value
      
    call FreeHost ( c_loc ( Value ) )
  
  end subroutine DeallocateHost_KDR_1D


  subroutine DeallocateHost_KDR_2D ( Value )
  
    real ( KDR ), dimension ( :, : ), target, intent ( in ) :: &
      Value
    
    call FreeHost ( c_loc ( Value ) )
    
  end subroutine DeallocateHost_KDR_2D


end module DeallocateHost_Command
