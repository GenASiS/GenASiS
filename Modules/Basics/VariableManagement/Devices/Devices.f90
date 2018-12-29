module Devices

  use AllocateDevice_Command
  use AssociateHost_Command
  use DeallocateDevice_Command
  use DisassociateHost_Command
  use UpdateDevice_Command
  use UpdateHost_Command
  use AllocateHost_Command
  use DeallocateHost_Command
  
  logical, parameter :: &
    DEVICE_TARGET = .false.

end module Devices
