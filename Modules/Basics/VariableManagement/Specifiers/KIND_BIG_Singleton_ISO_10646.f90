!--  KIND_BIG_Singleton defines kind parameters for bigger-than-normal 
!    intrinsic data types.

module KIND_BIG_Singleton

  use KIND_DEFAULT_Singleton

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      SUPERSCRIPT_1     = int ( z'00B9' ), &
      SUPERSCRIPT_2     = int ( z'00B2' ), &
      SUPERSCRIPT_3     = int ( z'00B3' ), &
      SUPERSCRIPT_4     = int ( z'2074' ), &
      SUPERSCRIPT_5     = int ( z'2075' ), &
      SUPERSCRIPT_6     = int ( z'2076' ), &
      SUPERSCRIPT_7     = int ( z'2077' ), &
      SUPERSCRIPT_8     = int ( z'2078' ), &
      SUPERSCRIPT_9     = int ( z'2079' ), &
      SUPERSCRIPT_MINUS = int ( z'207B' ), &
      CAPITAL_A_RING    = int ( z'00C5' ), &
      SMALL_H_STROKE    = int ( z'0127' ), &
      CIRCLE_DOT        = int ( z'2299' )

  type, public :: KindBigSingleton
    integer ( kind ( KDI ) ) :: &
      INTEGER           = selected_int_kind ( 15 ), &
      CHARACTER         = selected_char_kind ( 'ISO_10646' ), &
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
