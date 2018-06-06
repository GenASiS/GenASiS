module DisassociateDevice_Command

  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    DisassociateDevice
  
  interface DisassociateDevice
    module procedure DisassociateDevice_KDR_1D
    module procedure DisassociateDevice_KDR_2D
    module procedure DisassociateDevice_KDR_3D
  end interface DisassociateDevice
  
contains
  
  
  subroutine DisassociateDevice_KDR_1D ( Value, ErrorOption )
  
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      Error
    
    Error = DisassociateTarget ( c_loc ( Value ) )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine DisassociateDevice_KDR_1D


  subroutine DisassociateDevice_KDR_2D ( Value, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      Error
    
    Error = DisassociateTarget ( c_loc ( Value ) )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error

  end subroutine DisassociateDevice_KDR_2D
  
  
  subroutine DisassociateDevice_KDR_3D ( Value, ErrorOption )
  
    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      Error
    
    Error = DisassociateTarget ( c_loc ( Value ) )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error

  end subroutine DisassociateDevice_KDR_3D


end module DisassociateDevice_Command
