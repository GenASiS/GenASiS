program MeasuredValue_Form_Test

  use ISO_FORTRAN_ENV
  use KIND_DEFAULT_Singleton
  use KIND_BIG_Singleton
  use MeasuredValue_Form

  implicit none

  real ( KDR ) :: &
    A
  character ( 5 ) :: &
    Encoding
  type ( MeasuredValueForm ) :: &
    Length_1, &
    Length_2, &
    Time_1

!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  print *

  !-- Initializations

  call Length_1 % Initialize ( 'm', 10.0_KDR )
  call Length_2 % Initialize ( Length_1 ) 
  call Time_1 % Initialize ( 's', 0.5_KDR )

  print *
  print *, 'Length_1 =', Length_1
  print *, 'Length_2 =', Length_2
  print *, 'Time_1 =', Time_1

  !-- Addition

  print *
  print *, '*** Addition'
  print * 
  print *, 'Length_1 + Length_2 =', Length_1 + Length_2
  print *, 'Length_1 + Time_1 =', Length_1 + Time_1
  print *, 'Length_1 + 2.0 =', Length_1 + 2.0_KDR
  print *
  print *, '2.0 + Length_1 =', 2.0_KDR + Length_1
  print *
  print *, 'Length_1 + 3 =', Length_1 + 3_KDI
  print *
  print *, '3 + Length_1 =', 3_KDI + Length_1
  print *

  !-- Subtraction

  print *
  print *, '*** Subtraction'
  print *
  print *, 'Length_1 - Length_2 =', Length_1 - Length_2
  print *, 'Length_1 - Time_1 =', Length_1 - Time_1
  print *, 'Length_1 - 2.0 =', Length_1 - 2.0_KDR
  print *
  print *, '2.0 - Length_1 =', 2.0_KDR - Length_1
  print *
  print *, 'Length_1 - 3 =', Length_1 - 3_KDI
  print *
  print *, '3 - Length_1 =', 3_KDI - Length_1
  print *

  !-- Multiplication

  print *
  print *, '*** Multiplication'
  print *
  print *, 'Length_1 * Length_2 =', Length_1 * Length_2
  print *, 'Length_1 * Time_1 =', Length_1 * Time_1
  print *, 'Length_1 * 2.0 =', Length_1 * 2.0_KDR
  print *, '2.0 * Length_1 =', 2.0_KDR * Length_1
  print *, 'Length_1 * 3 =', Length_1 * 3_KDI
  print *, '3 * Length_1 =', 3_KDI * Length_1

  !-- Exponentiation

  print *
  print *, '*** Exponentiation'
  print *
  print *, 'Length_1 ** 2 =', Length_1 ** 2_KDI
  print *, 'Length_1 ** (-3) =', Length_1 ** (-3_KDI)
  print *, '( Length_1 ** 2 ) ** (-3)  =', &
           ( Length_1 ** 2_KDI ) ** (-3_KDI)
!  print *, 'Length_1 ** 0.5 =', Length_1 ** 0.5_KDR

  !-- Division

  print *
  print *, '*** Division'
  print *
  print *, 'Length_1 / Length_2 =', Length_1 / Length_2
  print *, 'Length_1 / Time_1 =', Length_1 / Time_1
  print *, 'Length_1 / 2.0 =', Length_1 / 2.0_KDR
  print *, '2.0 / Length_1 =', 2.0_KDR / Length_1
  print *, 'Length_1 / 3 =', Length_1 / 3_KDI
  print *, '3 / Length_1 =', 3_KDI / Length_1

  !-- Assignment

  print *
  print *, '*** Assignment'
  print *
  A = Length_1
  print *, 'A = Length_1 =', A
  print *

  !-- EqualTo

  print *
  print *, '*** EqualTo'
  print *
  print *, 'Length_1  ==  Length_2 =', Length_1 == Length_2
  print *
  print *, '2 * Length_1  ==  Length_2 =', 2 * Length_1 == Length_2
  print *
  print *, 'Length_1  ==  2 * Length_2 =', Length_1 == 2 * Length_2
  print *
  print *, 'Length_1  ==  1.0 =', Length_1 == 1.0_KDR
  print *
  print *, 'Length_1  ==  10.0 =', Length_1 == 10.0_KDR
  print *
  print *, 'Length_1  ==  100.0 =', Length_1 == 100.0_KDR
  print *
  print *, '1.0  ==  Length_1 =', 1.0_KDR == Length_1
  print *
  print *, '10.0  ==  Length_1 =', 10.0_KDR == Length_1
  print *
  print *, '100.0  ==  Length_1 =', 100.0_KDR == Length_1
  print *
  print *, 'Length_1  ==  1 =', Length_1 == 1_KDI
  print *
  print *, 'Length_1  ==  10 =', Length_1 == 10_KDI
  print *
  print *, 'Length_1  ==  100 =', Length_1 == 100_KDI
  print *
  print *, '1  ==  Length_1 =', 1_KDI == Length_1
  print *
  print *, '10  ==  Length_1 =', 10_KDI == Length_1
  print *
  print *, '100  ==  Length_1 =', 100_KDI == Length_1
  print *

  !-- NotEqualTo

  print *
  print *, '*** NotEqualTo'
  print *
  print *, 'Length_1  /=  Length_2 =', Length_1 /= Length_2
  print *
  print *, '2 * Length_1  /=  Length_2 =', 2 * Length_1 /= Length_2
  print *
  print *, 'Length_1  /=  2 * Length_2 =', Length_1 /= 2 * Length_2
  print *
  print *, 'Length_1  /=  1.0 =', Length_1 /= 1.0_KDR
  print *
  print *, 'Length_1  /=  10.0 =', Length_1 /= 10.0_KDR
  print *
  print *, 'Length_1  /=  100.0 =', Length_1 /= 100.0_KDR
  print *
  print *, '1.0  /=  Length_1 =', 1.0_KDR /= Length_1
  print *
  print *, '10.0  /=  Length_1 =', 10.0_KDR /= Length_1
  print *
  print *, '100.0  /=  Length_1 =', 100.0_KDR /= Length_1
  print *
  print *, 'Length_1  /=  1 =', Length_1 /= 1_KDI
  print *
  print *, 'Length_1  /=  10 =', Length_1 /= 10_KDI
  print *
  print *, 'Length_1  /=  100 =', Length_1 /= 100_KDI
  print *
  print *, '1  /=  Length_1 =', 1_KDI /= Length_1
  print *
  print *, '10  /=  Length_1 =', 10_KDI /= Length_1
  print *
  print *, '100  /=  Length_1 =', 100_KDI /= Length_1
  print *

  !-- GreaterThan

  print *
  print *, '*** GreaterThan'
  print *
  print *, 'Length_1  >  Length_2 =', Length_1 > Length_2
  print *
  print *, '2 * Length_1  >  Length_2 =', 2 * Length_1 > Length_2
  print *
  print *, 'Length_1  >  2 * Length_2 =', Length_1 > 2 * Length_2
  print *
  print *, 'Length_1  >  1.0 =', Length_1 > 1.0_KDR
  print *
  print *, 'Length_1  >  10.0 =', Length_1 > 10.0_KDR
  print *
  print *, 'Length_1  >  100.0 =', Length_1 > 100.0_KDR
  print *
  print *, '1.0  >  Length_1 =', 1.0_KDR > Length_1
  print *
  print *, '10.0  >  Length_1 =', 10.0_KDR > Length_1
  print *
  print *, '100.0  >  Length_1 =', 100.0_KDR > Length_1
  print *
  print *, 'Length_1  >  1 =', Length_1 > 1_KDI
  print *
  print *, 'Length_1  >  10 =', Length_1 > 10_KDI
  print *
  print *, 'Length_1  >  100 =', Length_1 > 100_KDI
  print *
  print *, '1  >  Length_1 =', 1_KDI > Length_1
  print *
  print *, '10  >  Length_1 =', 10_KDI > Length_1
  print *
  print *, '100  >  Length_1 =', 100_KDI > Length_1
  print *

  !-- LessThan

  print *
  print *, '*** LessThan'
  print *
  print *, 'Length_1  <  Length_2 =', Length_1 < Length_2
  print *
  print *, '2 * Length_1  <  Length_2 =', 2 * Length_1 < Length_2
  print *
  print *, 'Length_1  <  2 * Length_2 =', Length_1 < 2 * Length_2
  print *
  print *, 'Length_1  <  1.0 =', Length_1 < 1.0_KDR
  print *
  print *, 'Length_1  <  10.0 =', Length_1 < 10.0_KDR
  print *
  print *, 'Length_1  <  100.0 =', Length_1 < 100.0_KDR
  print *
  print *, '1.0  <  Length_1 =', 1.0_KDR < Length_1
  print *
  print *, '10.0  <  Length_1 =', 10.0_KDR < Length_1
  print *
  print *, '100.0  <  Length_1 =', 100.0_KDR < Length_1
  print *
  print *, 'Length_1  <  1 =', Length_1 < 1_KDI
  print *
  print *, 'Length_1  <  10 =', Length_1 < 10_KDI
  print *
  print *, 'Length_1  <  100 =', Length_1 < 100_KDI
  print *
  print *, '1  <  Length_1 =', 1_KDI < Length_1
  print *
  print *, '10  <  Length_1 =', 10_KDI < Length_1
  print *
  print *, '100  <  Length_1 =', 100_KDI < Length_1
  print *

  !-- GreaterThanEqualTo

  print *
  print *, '*** GreaterThanEqualTo'
  print *
  print *, 'Length_1  >=  Length_2 =', Length_1 >= Length_2
  print *
  print *, '2 * Length_1  >=  Length_2 =', 2 * Length_1 >= Length_2
  print *
  print *, 'Length_1  >=  2 * Length_2 =', Length_1 >= 2 * Length_2
  print *
  print *, 'Length_1  >=  1.0 =', Length_1 >= 1.0_KDR
  print *
  print *, 'Length_1  >=  10.0 =', Length_1 >= 10.0_KDR
  print *
  print *, 'Length_1  >=  100.0 =', Length_1 >= 100.0_KDR
  print *
  print *, '1.0  >=  Length_1 =', 1.0_KDR >= Length_1
  print *
  print *, '10.0  >=  Length_1 =', 10.0_KDR >= Length_1
  print *
  print *, '100.0  >=  Length_1 =', 100.0_KDR >= Length_1
  print *
  print *, 'Length_1  >=  1 =', Length_1 >= 1_KDI
  print *
  print *, 'Length_1  >=  10 =', Length_1 >= 10_KDI
  print *
  print *, 'Length_1  >=  100 =', Length_1 >= 100_KDI
  print *
  print *, '1  >=  Length_1 =', 1_KDI >= Length_1
  print *
  print *, '10  >=  Length_1 =', 10_KDI >= Length_1
  print *
  print *, '100  >=  Length_1 =', 100_KDI >= Length_1
  print *

  !-- LessThanEqualTo

  print *
  print *, '*** LessThanEqualTo'
  print *
  print *, 'Length_1  <=  Length_2 =', Length_1 <= Length_2
  print *
  print *, '2 * Length_1  <=  Length_2 =', 2 * Length_1 <= Length_2
  print *
  print *, 'Length_1  <=  2 * Length_2 =', Length_1 <= 2 * Length_2
  print *
  print *, 'Length_1  <=  1.0 =', Length_1 <= 1.0_KDR
  print *
  print *, 'Length_1  <=  10.0 =', Length_1 <= 10.0_KDR
  print *
  print *, 'Length_1  <=  100.0 =', Length_1 <= 100.0_KDR
  print *
  print *, '1.0  <=  Length_1 =', 1.0_KDR <= Length_1
  print *
  print *, '10.0  <=  Length_1 =', 10.0_KDR <= Length_1
  print *
  print *, '100.0  <=  Length_1 =', 100.0_KDR <= Length_1
  print *
  print *, 'Length_1  <=  1 =', Length_1 <= 1_KDI
  print *
  print *, 'Length_1  <=  10 =', Length_1 <= 10_KDI
  print *
  print *, 'Length_1  <=  100 =', Length_1 <= 100_KDI
  print *
  print *, '1  <=  Length_1 =', 1_KDI <= Length_1
  print *
  print *, '10  <=  Length_1 =', 10_KDI <= Length_1
  print *
  print *, '100  <=  Length_1 =', 100_KDI <= Length_1
  print *

end program MeasuredValue_Form_Test
