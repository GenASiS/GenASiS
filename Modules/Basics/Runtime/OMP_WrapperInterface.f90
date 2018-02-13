module OMP_WrapperInterface

  use ISO_C_BINDING
  
  implicit none
  private
  
  public :: &
    Allocate_D, &
    AssociateTarget, &
    DisassociateTarget
  
  interface 

    type ( c_ptr ) function Allocate_D ( nValues ) &
                              bind ( c, name = 'Allocate_D_Double' )
      use iso_c_binding
      implicit none
      integer ( kind = c_int ), value :: &
        nValues
    end function Allocate_D
    
    
    integer ( c_int ) function AssociateTarget &
                     ( Host, Device, nValues, oValue ) &
                     bind ( c, name = 'AssociateTargetWrapper' )
    
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host, &
        Device
      integer ( c_int ), value :: &
        nValues, &
        oValue
    end function AssociateTarget

    
    integer ( c_int ) function DisassociateTarget ( Host ) &
                     bind ( c, name = 'DisassociateTargetWrapper' )
    
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Host
    end function DisassociateTarget

  end interface 

end module OMP_WrapperInterface
