module Device_C

  use iso_c_binding
  
  implicit none
  private
  
  public :: &
    OnTarget, &
    AllocateTargetDouble, &
    AssociateTargetDouble, &
    DeallocateTarget, &
    DisassociateTarget, &
    AllocateHostDouble, &
    FreeHost, &
    HostToDeviceCopyDouble, &
    DeviceToHostCopyDouble, &
    IsOffloadEnabled
  
  interface 
  
    integer ( c_int ) function OnTarget ( Host ) &
                                 bind ( c, name = 'OnTarget_OMP' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host
    end function OnTarget


    type ( c_ptr ) function AllocateTargetDouble ( nValues ) &
                              bind ( c, name = 'AllocateTargetDouble_OMP' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        nValues
    end function AllocateTargetDouble
    
    
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
                              bind ( c, name = 'AllocateHostDouble_CUDA' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        nValues
    end function AllocateHostDouble
    
    
    subroutine FreeHost ( Host ) bind ( c, name = 'FreeHost_CUDA' )
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

  end interface 

end module Device_C
