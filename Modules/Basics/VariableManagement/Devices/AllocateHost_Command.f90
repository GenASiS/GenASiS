module AllocateHost_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    AllocateHost
  
  interface AllocateHost
    module procedure AllocateHost_KDR_1D
    module procedure AllocateHost_KDR_2D
  end interface AllocateHost
  
contains

  
  subroutine AllocateHost_KDR_1D ( Value, Shape )
  
    real ( KDR ), dimension ( : ), pointer, intent ( inout ) :: &
      Value
    integer ( KDI ), dimension ( 1 ), intent ( in ) :: &
      Shape
      
    type ( c_ptr ) :: &
      Address
    
    Address = AllocateHostDouble ( Shape ( 1 ) )
    call c_f_pointer ( Address, Value, Shape )
  
  end subroutine AllocateHost_KDR_1D


  subroutine AllocateHost_KDR_2D ( Value, Shape )
  
    real ( KDR ), dimension ( :, : ), pointer, intent ( inout ) :: &
      Value
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      Shape
    
    type ( c_ptr ) :: &
      Address
    
    Address = AllocateHostDouble ( product ( Shape ) )
    call c_f_pointer ( Address, Value, Shape )
  
  end subroutine AllocateHost_KDR_2D


end module AllocateHost_Command
