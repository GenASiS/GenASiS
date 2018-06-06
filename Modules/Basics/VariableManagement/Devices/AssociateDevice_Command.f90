module AssociateDevice_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    AssociateDevice
  
  interface AssociateDevice
    module procedure AssociateDevice_KDR_1D
    module procedure AssociateDevice_KDR_2D
    module procedure AssociateDevice_KDR_3D
  end interface AssociateDevice
  
contains
  
  
  subroutine AssociateDevice_KDR_1D &
               ( Device, Value, oValueOption, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      oValueOption
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      oValue, &
      Error
    
    oValue = 0
    if ( present ( oValueOption ) ) &
      oValue = oValueOption 
          
    Error = AssociateTargetDouble &
              ( c_loc ( Value ), Device, size ( Value ), oValue )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine AssociateDevice_KDR_1D


  subroutine AssociateDevice_KDR_2D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      oValue, &
      Error
    
    oValue = 0
          
    Error = AssociateTargetDouble &
              ( c_loc ( Value ), Device, size ( Value ), oValue )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error

  end subroutine AssociateDevice_KDR_2D

  
  subroutine AssociateDevice_KDR_3D ( Device, Value, ErrorOption )
  
    type ( c_ptr ), intent ( in ) :: &
      Device
    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      oValue, &
      Error
    
    oValue = 0
          
    Error = AssociateTargetDouble &
              ( c_loc ( Value ), Device, size ( Value ), oValue )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error

  end subroutine AssociateDevice_KDR_3D


end module AssociateDevice_Command
