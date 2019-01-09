!--  KIND_TINY_Singleton defines kind parameters for tinier-than-normal 
!    intrinsic data types.

module KIND_TINY_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

  type, public :: KindTinySingleton
    integer ( kind ( KDI ) ) :: &
      INTEGER = selected_int_kind ( 2 ), &
      LOGICAL = selected_int_kind ( 1 )  !-- No selected_logical_kind!
  end type KindTinySingleton

  type ( KindTinySingleton ), public, parameter :: &
    KIND_TINY = KindTinySingleton ( )

  integer ( KDI ), public, parameter :: &
    KTI = KIND_TINY % INTEGER, &
    KTL = KIND_TINY % LOGICAL

end module KIND_TINY_Singleton
