!--  KIND_BIG_Singleton defines kind parameters for bigger-than-normal 
!    intrinsic data types.

module KIND_BIG_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

    !-- For ASCII, replacing unicode characters with "?"; in principle
    !   they should never be seen.
    integer ( KDI ), private, parameter :: &
      SUPERSCRIPT_1     = int ( z'003F' ), &
      SUPERSCRIPT_2     = int ( z'003F' ), &
      SUPERSCRIPT_3     = int ( z'003F' ), &
      SUPERSCRIPT_4     = int ( z'003F' ), &
      SUPERSCRIPT_5     = int ( z'003F' ), &
      SUPERSCRIPT_6     = int ( z'003F' ), &
      SUPERSCRIPT_7     = int ( z'003F' ), &
      SUPERSCRIPT_8     = int ( z'003F' ), &
      SUPERSCRIPT_9     = int ( z'003F' ), &
      SUPERSCRIPT_MINUS = int ( z'003F' ), &
      CAPITAL_A_RING    = int ( z'003F' ), &
      SMALL_H_STROKE    = int ( z'003F' ), &
      CIRCLE_DOT        = int ( z'003F' )

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
      SUPERSCRIPT_MINUS = SUPERSCRIPT_MINUS, &
      CAPITAL_A_RING    = CAPITAL_A_RING, &
      SMALL_H_STROKE    = SMALL_H_STROKE, &
      CIRCLE_DOT        = CIRCLE_DOT
  end type KindBigSingleton

  type ( KindBigSingleton ), public, parameter :: &
    KIND_BIG = KindBigSingleton ( )

  integer ( KDI ), public, parameter :: &
    KBI  = KIND_BIG % INTEGER, &
    KBCH = KIND_BIG % CHARACTER

end module KIND_BIG_Singleton
