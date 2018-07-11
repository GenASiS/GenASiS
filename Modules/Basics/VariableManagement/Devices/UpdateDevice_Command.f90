module UpdateDevice_Command
  
  use iso_c_binding
  use Specifiers
  use AssociateDevice_Command
  use DisassociateDevice_Command
  
  implicit none
  private
  
  public :: &
    UpdateDevice, &
    UpdateDeviceAsync, &
    FinishUpdate
  
  interface UpdateDevice
    module procedure UpdateDevice_KDR_1D
    module procedure UpdateDevice_KDR_2D
  end interface UpdateDevice
  
  interface UpdateDeviceAsync
    module procedure UpdateDevice_KDR_1D_Async
    module procedure UpdateDevice_KDR_2D_Async
  end interface UpdateDeviceAsync
  
contains

  
  subroutine UpdateDevice_KDR_1D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateDevice ( Device, Value, ErrorOption )
    !$OMP target update to ( Value ) 
    call DisassociateDevice ( Value, ErrorOption )
  
  end subroutine UpdateDevice_KDR_1D


  subroutine UpdateDevice_KDR_2D ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateDevice ( Device, Value, ErrorOption )
    !$OMP target update to ( Value ) 
    call DisassociateDevice ( Value, ErrorOption )
  
  end subroutine UpdateDevice_KDR_2D


  subroutine UpdateDevice_KDR_1D_Async ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateDevice ( Device, Value, ErrorOption )
    !$OMP target update to ( Value ) nowait depend ( out: Device )
    call DisassociateDevice ( Value, ErrorOption )
  
  end subroutine UpdateDevice_KDR_1D_Async


  subroutine UpdateDevice_KDR_2D_Async ( Value, Device, ErrorOption )
  
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    type ( c_ptr ), intent ( in ) :: &
      Device
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
      
    call AssociateDevice ( Device, Value, ErrorOption )
    !$OMP target update to ( Value ) nowait depend ( out: Device )
    call DisassociateDevice ( Value, ErrorOption )
  
  end subroutine UpdateDevice_KDR_2D_Async
  
  
  subroutine FinishUpdate ( Device )
    
    type ( c_ptr ), intent ( in ) :: &
      Device
      
    integer ( KDI ) :: &
      Dummy
    
    !$OMP task depend ( in: Device )
    Dummy = 1
    !$OMP end task
  
  
  end subroutine FinishUpdate


end module UpdateDevice_Command
