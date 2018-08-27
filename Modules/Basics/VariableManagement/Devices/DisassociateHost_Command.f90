module DisassociateHost_Command

  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    DisassociateHost
  
  interface DisassociateHost
    module procedure DisassociateHost_KDR_1D
    module procedure DisassociateHost_KDR_2D
    module procedure DisassociateHost_KDR_3D
  end interface DisassociateHost
  
contains
  
  
  subroutine DisassociateHost_KDR_1D ( Value, ErrorOption )
  
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      Error
    
    Error = DisassociateTarget ( c_loc ( Value ) )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
  
  end subroutine DisassociateHost_KDR_1D


  subroutine DisassociateHost_KDR_2D ( Value, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      Error
    
    Error = DisassociateTarget ( c_loc ( Value ) )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error

  end subroutine DisassociateHost_KDR_2D
  
  
  subroutine DisassociateHost_KDR_3D ( Value, ErrorOption )
  
    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
      Value
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    
    integer ( c_int ) :: &
      Error
    
    Error = DisassociateTarget ( c_loc ( Value ) )
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error

  end subroutine DisassociateHost_KDR_3D


end module DisassociateHost_Command
