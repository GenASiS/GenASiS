!--  KIND_BIG_Singleton defines kind parameters for bigger-than-normal 
!    intrinsic data types.

module KIND_BIG_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

  type, public :: KindBigSingleton
    integer ( kind ( KDI ) ) :: &
      INTEGER = selected_int_kind ( 15 )
  end type KindBigSingleton

  type ( KindBigSingleton ), public, parameter :: &
    KIND_BIG = KindBigSingleton ( )

  integer ( KDI ), public, parameter :: &
    KBI = KIND_BIG % INTEGER

end module KIND_BIG_Singleton
