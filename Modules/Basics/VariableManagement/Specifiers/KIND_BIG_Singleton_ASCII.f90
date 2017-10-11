!--  KIND_BIG_Singleton defines kind parameters for bigger-than-normal 
!    intrinsic data types.

module KIND_BIG_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

    !-- Safe values for ASCII; normal digits
    integer ( KDI ), private, parameter :: &
      SUPERSCRIPT_1     = int ( z'0031' ), &
      SUPERSCRIPT_2     = int ( z'0032' ), &
      SUPERSCRIPT_3     = int ( z'0033' ), &
      SUPERSCRIPT_4     = int ( z'0034' ), &
      SUPERSCRIPT_5     = int ( z'0035' ), &
      SUPERSCRIPT_6     = int ( z'0036' ), &
      SUPERSCRIPT_7     = int ( z'0037' ), &
      SUPERSCRIPT_8     = int ( z'0038' ), &
      SUPERSCRIPT_9     = int ( z'0039' ), &
      SUPERSCRIPT_MINUS = int ( z'002D' )

  type, public :: KindBigSingleton
    integer ( kind ( KDI ) ) :: &
      INTEGER           = selected_int_kind ( 15 ), &
      CHARACTER         = selected_char_kind ( 'ASCII' ), &
      SUPERSCRIPT_1     = SUPERSCRIPT_1, &
      SUPERSCRIPT_2     = SUPERSCRIPT_2, &
      SUPERSCRIPT_3     = SUPERSCRIPT_3, &
      SUPERSCRIPT_4     = SUPERSCRIPT_4, &
      SUPERSCRIPT_5     = SUPERSCRIPT_5, &
      SUPERSCRIPT_6     = SUPERSCRIPT_6, &
      SUPERSCRIPT_7     = SUPERSCRIPT_7, &
      SUPERSCRIPT_8     = SUPERSCRIPT_8, &
      SUPERSCRIPT_9     = SUPERSCRIPT_9, &
      SUPERSCRIPT_MINUS = SUPERSCRIPT_MINUS
  end type KindBigSingleton

  type ( KindBigSingleton ), public, parameter :: &
    KIND_BIG = KindBigSingleton ( )

  integer ( KDI ), public, parameter :: &
    KBI  = KIND_BIG % INTEGER, &
    KBCH = KIND_BIG % CHARACTER

end module KIND_BIG_Singleton
