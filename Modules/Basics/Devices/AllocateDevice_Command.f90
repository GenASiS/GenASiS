module AllocateDevice_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    AllocateDevice
  
  interface AllocateDevice
    module procedure AllocateDevice_KDI_1D
    module procedure AllocateDevice_KDR
    module procedure AllocateDevice_KDR_1D
    module procedure AllocateDevice_KDR_2D
    module procedure AllocateDevice_KDR_4D
  end interface AllocateDevice
  
contains


  subroutine AllocateDevice_KDI_1D ( Value, Device )
  
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( out ) :: &
      Device
    
    Device = AllocateTargetInteger ( size ( Value ) )
  
  end subroutine AllocateDevice_KDI_1D
  
  
  subroutine AllocateDevice_KDR ( nValues, Device )
  
    integer ( KDI ), intent ( in ) :: &
      nValues
    type ( c_ptr ), intent ( out ) :: &
      Device
    
    Device = AllocateTargetDouble ( nValues )
  
  end subroutine AllocateDevice_KDR


  subroutine AllocateDevice_KDR_1D ( Value, Device )
  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( out ) :: &
      Device
      
    Device = AllocateTargetDouble ( size ( Value ) )
  
  end subroutine AllocateDevice_KDR_1D


  subroutine AllocateDevice_KDR_2D ( Value, Device )
  
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( out ) :: &
      Device
    
    Device = AllocateTargetDouble ( size ( Value ) )
  
  end subroutine AllocateDevice_KDR_2D


  subroutine AllocateDevice_KDR_4D ( Value, Device )
  
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( out ) :: &
      Device
    
    Device = AllocateTargetDouble ( size ( Value ) )
  
  end subroutine AllocateDevice_KDR_4D


end module AllocateDevice_Command
