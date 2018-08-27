module AssociateHost_Command
  
  use iso_c_binding
  use Specifiers
  use Device_C
  
  implicit none
  private
  
  public :: &
    AssociateHost
  
  interface AssociateHost
    module procedure AssociateHost_KDR_1D
    module procedure AssociateHost_KDR_2D
    module procedure AssociateHost_KDR_3D
  end interface AssociateHost
  
contains
  
  
  subroutine AssociateHost_KDR_1D &
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
  
  end subroutine AssociateHost_KDR_1D


  subroutine AssociateHost_KDR_2D ( Device, Value, ErrorOption )
  
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

  end subroutine AssociateHost_KDR_2D

  
  subroutine AssociateHost_KDR_3D ( Device, Value, ErrorOption )
  
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

  end subroutine AssociateHost_KDR_3D


end module AssociateHost_Command
