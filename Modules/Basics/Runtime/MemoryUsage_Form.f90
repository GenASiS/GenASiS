module MemoryUsage_Form
  
  use Specifiers
  use Devices
  use Display
  use MessagePassing
  use GetMemoryUsage_Command

  implicit none
  private
  
  type, public :: MemoryUsageForm
    type ( QuantityForm )  :: &
      HighWaterMark, &
      HighWaterMarkMax, &
      HighWaterMarkMin, &
      HighWaterMarkMean
    type ( QuantityForm )  :: &
      ResidentSetSize, &
      ResidentSetSizeMax, &
      ResidentSetSizeMin, &
      ResidentSetSizeMean
    type ( QuantityForm )  :: &
      DeviceMemoryTotal, &
      DeviceMemoryFree, &
      DeviceMemoryUsed
  contains
    procedure, public, pass :: &
      Compute
  end type MemoryUsageForm


contains


  subroutine Compute ( MU, Ignorability, CommunicatorOption )

    class ( MemoryUsageForm ), intent ( inout ) :: &
      MU
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption

    if ( present ( CommunicatorOption ) ) then
      call GetMemoryUsage &
             ( MU % HighWaterMark, MU % ResidentSetSize, Ignorability, &
               C_Option        = CommunicatorOption, &
               Max_HWM_Option  = MU % HighWaterMarkMax, &
               Min_HWM_Option  = MU % HighWaterMarkMin, &
               Mean_HWM_Option = MU % HighWaterMarkMean, &
               Max_RSS_Option  = MU % ResidentSetSizeMax, &
               Min_RSS_Option  = MU % ResidentSetSizeMin, &
               Mean_RSS_Option = MU % ResidentSetSizeMean )
    else
      call GetMemoryUsage &
             ( MU % HighWaterMark, MU % ResidentSetSize, Ignorability )
    end if
    
    call Show ( MU % HighWaterMark, 'This process HWM', Ignorability )
    call Show ( MU % ResidentSetSize, 'This process RSS', Ignorability )
    
    if ( present ( CommunicatorOption ) ) then
      
      call Show ( MU % HighWaterMarkMax, 'Across processes max HWM', &
                  Ignorability )
      call Show ( MU % HighWaterMarkMin, 'Across processes min HWM', &
                  Ignorability )
      call Show ( MU % HighWaterMarkMean, 'Across processes mean HWM', &
                  Ignorability )
      
      call Show ( MU % ResidentSetSizeMax, 'Across processes max RSS', &
                  Ignorability )
      call Show ( MU % ResidentSetSizeMin, 'Across processes min RSS', &
                  Ignorability )
      call Show ( MU % ResidentSetSizeMean, 'Across processes mean RSS', &
                  Ignorability )
    
    end if
    
    if ( OffloadEnabled ( ) ) then
      call Show ( 'Device memory info', Ignorability )
      call GetDeviceMemoryInfo &
             ( MU % DeviceMemoryTotal, MU % DeviceMemoryUsed, &
               MU % DeviceMemoryFree )
      call Show ( MU % DeviceMemoryTotal, 'This process device memory', &
                  Ignorability )
      call Show ( MU % DeviceMemoryUsed,  'This process device memory used', &
                  Ignorability )
      call Show ( MU % DeviceMemoryFree,  'This process device memory free', &
                  Ignorability )
    end if
    
  end subroutine Compute


end module MemoryUsage_Form
