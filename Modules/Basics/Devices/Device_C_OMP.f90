module Device_C

  use iso_c_binding
  
  implicit none
  private
  
  public :: &
    SetDevice, &
    GetDevice, &
    OnTarget, &
    AllocateTargetInteger, &
    AllocateTargetDouble, &
    AllocateTargetLogical, &
    AssociateTargetInteger, &
    AssociateTargetDouble, &
    AssociateTargetLogical, &
    DeallocateTarget, &
    DisassociateTarget, &
    AllocateHostDouble, &
    FreeHost, &
    HostToDeviceCopyDouble, &
    DeviceToHostCopyDouble, &
    IsOffloadEnabled, &
    DeviceMemGetInfo
  
  interface 
    
    integer ( c_int ) function SetDevice ( iDevice ) &
                                 bind ( c, name = 'SetDevice' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        iDevice
    end function SetDevice
    

    integer ( c_int ) function GetDevice ( iDevice ) &
                                 bind ( c, name = 'GetDevice' )
      use iso_c_binding
      implicit none
      integer ( c_int ) :: &
        iDevice
    end function GetDevice
    
    
    integer ( c_int ) function OnTarget ( Host ) &
                                 bind ( c, name = 'OnTarget_OMP' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host
    end function OnTarget


    type ( c_ptr ) function AllocateTargetInteger ( nValues ) &
                              bind ( c, name = 'AllocateTargetInteger_OMP' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        nValues
    end function AllocateTargetInteger
    
    
    type ( c_ptr ) function AllocateTargetDouble ( nValues ) &
                              bind ( c, name = 'AllocateTargetDouble_OMP' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        nValues
    end function AllocateTargetDouble
    
    
    type ( c_ptr ) function AllocateTargetLogical ( nValues ) &
                              bind ( c, name = 'AllocateTargetLogical_OMP' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        nValues
    end function AllocateTargetLogical
    
    
    integer ( c_int ) function AssociateTargetInteger &
                        ( Host, Device, nValues, oValue ) &
                        bind ( c, name = 'AssociateTargetInteger_OMP' )
    
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host, &
        Device
      integer ( c_int ), value :: &
        nValues, &
        oValue
    end function AssociateTargetInteger
    
    
    integer ( c_int ) function AssociateTargetDouble &
                        ( Host, Device, nValues, oValue ) &
                        bind ( c, name = 'AssociateTargetDouble_OMP' )
    
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host, &
        Device
      integer ( c_int ), value :: &
        nValues, &
        oValue
    end function AssociateTargetDouble
    
    
    integer ( c_int ) function AssociateTargetLogical &
                        ( Host, Device, nValues, oValue ) &
                        bind ( c, name = 'AssociateTargetLogical_OMP' )
    
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host, &
        Device
      integer ( c_int ), value :: &
        nValues, &
        oValue
    end function AssociateTargetLogical
    
    
    subroutine DeallocateTarget ( Device ) bind ( c, name = 'FreeTarget_OMP' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Device
    
    end subroutine DeallocateTarget

    
    integer ( c_int ) function DisassociateTarget ( Host ) &
                        bind ( c, name = 'DisassociateTarget_OMP' )
    
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host
    end function DisassociateTarget
    
    
    type ( c_ptr ) function AllocateHostDouble ( nValues ) &
                              bind ( c, name = 'AllocateHostDouble_Device' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        nValues
    end function AllocateHostDouble
    
    
    subroutine FreeHost ( Host ) bind ( c, name = 'FreeHost_Device' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host
    
    end subroutine FreeHost
    
    
    integer ( c_int ) function HostToDeviceCopyDouble &
                        ( Host, Device, nValues, oHostValue, &
                          oDeviceValue ) &
                          bind ( c, name = 'HostToDeviceCopyDouble_OMP' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host, &
        Device
      integer ( c_int ), value :: &
        nValues, &
        oHostValue, &
        oDeviceValue
    end function HostToDeviceCopyDouble

    
    integer ( c_int ) function DeviceToHostCopyDouble &
                        ( Device, Host, nValues, oDeviceValue, &
                          oHostValue ) &
                          bind ( c, name = 'DeviceToHostCopyDouble_OMP' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Device, &
        Host
      integer ( c_int ), value :: &
        nValues, &
        oDeviceValue, &
        oHostValue
    end function DeviceToHostCopyDouble
    
    
    logical ( c_bool ) function IsOffloadEnabled ( ) &
                           bind ( c, name = 'OffloadEnabled' )
      use iso_c_binding
      implicit none
    end function IsOffloadEnabled
    
    
    integer ( c_int ) function DeviceMemGetInfo ( Free, Total ) &
                        bind ( c, name = 'DeviceMemGetInfo_Device' )
      use iso_c_binding
      implicit none
      integer ( c_size_t ) :: &
        Free, &
        Total
    end function DeviceMemGetInfo 
    
  end interface 

end module Device_C
