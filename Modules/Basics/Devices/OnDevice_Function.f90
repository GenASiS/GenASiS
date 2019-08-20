module OnDevice_Function

  use iso_c_binding
  use Specifiers
  use Device_C

  implicit none
  private

  interface OnDevice
    module procedure OnDevice_1D
    module procedure OnDevice_3D
  end interface OnDevice

  public :: &
    OnDevice

contains


  function OnDevice_1D ( Value ) result ( IP )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Value

    integer ( c_int ) :: &
      Retval
    logical :: &
      IP

    Retval = OnTarget ( c_loc ( Value ) )

    if ( Retval /= 0 ) then
      IP = .true.
    else
      IP = .false.
    end if

  end function OnDevice_1D


  function OnDevice_3D ( Value ) result ( IP )

    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
      Value

    integer ( c_int ) :: &
      Retval
    logical :: &
      IP

    Retval = OnTarget ( c_loc ( Value ) )

    if ( Retval /= 0 ) then
      IP = .true.
    else
      IP = .false.
    end if

  end function OnDevice_3D


end module OnDevice_Function
