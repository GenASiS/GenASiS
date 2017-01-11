!-- KIND_DEFAULT_Singleton defines default kind parameters for
!   intrinsic data types

module KIND_DEFAULT_Singleton

  implicit none
  private

  type, public :: KindDefaultSingleton
    integer ( kind ( 1 ) ) :: &
      INTEGER = kind ( 1 ), &
      REAL    = kind ( 1.d0 ), &
      COMPLEX = kind ( ( 1.d0, 1.d0 ) ), &
      LOGICAL = kind ( .true. )
  end type KindDefaultSingleton

  type ( KindDefaultSingleton ), public, parameter :: &
    KIND_DEFAULT = KindDefaultSingleton ( )

  integer ( KIND_DEFAULT % INTEGER ), public, parameter :: &
    KDI = KIND_DEFAULT % INTEGER, &
    KDR = KIND_DEFAULT % REAL, &
    KDC = KIND_DEFAULT % COMPLEX, &
    KDL = KIND_DEFAULT % LOGICAL

end module KIND_DEFAULT_Singleton
