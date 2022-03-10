module MemoryUsage_C_macOS

  use iso_c_binding
  
  implicit none
  private

  public :: &
    ShowMemory_macOS, &
    GetMemory_macOS

  interface

    subroutine ShowMemory_macOS ( ) bind ( c, name = 'memcheck' )
    end subroutine ShowMemory_macOS

    subroutine GetMemory_macOS ( RSS_kB_cloc ) bind ( c, name = 'memget' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        RSS_kB_cloc
    end subroutine GetMemory_macOS

  end interface

end module MemoryUsage_C_macOS
