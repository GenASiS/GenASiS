!--  KIND_SMALL_Singleton defines kind parameters for smaller-than-normal 
!    intrinsic data types.

module KIND_SMALL_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

  type, public :: KindSmallSingleton
    integer ( kind ( KDI ) ) :: &
      INTEGER = selected_int_kind ( 4 )
  end type KindSmallSingleton

  type ( KindSmallSingleton ), public, parameter :: &
    KIND_SMALL = KindSmallSingleton ( )

  integer ( KDI ), public, parameter :: &
    KSI = KIND_SMALL % INTEGER

end module KIND_SMALL_Singleton
