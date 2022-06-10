module GetMemoryUsage_macOS_Command

  use iso_c_binding
  
  implicit none
  private

  public :: &
    GetMemoryUsage_macOS

  interface

    subroutine GetMemoryUsage_macOS ( HWM_kB_cloc, RSS_kB_cloc ) &
      bind ( c, name = 'GetMemoryUsage_macOS_C' )
        use iso_c_binding
        implicit none
        type ( c_ptr ), value :: &
          HWM_kB_cloc, &
          RSS_kB_cloc
    end subroutine GetMemoryUsage_macOS

  end interface

end module GetMemoryUsage_macOS_Command
