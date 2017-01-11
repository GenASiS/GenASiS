program ConsoleHeader_Form_Test

  use ConsoleHeader_Form

  implicit none

  type ( ConsoleHeaderForm ) :: &
    CH

  print *, trim ( CH % LABEL ( 1 ) ) // ' = ', CH % ERROR
  print *, trim ( CH % LABEL ( 2 ) ) // ' = ', CH % WARNING
  print *, trim ( CH % LABEL ( 3 ) ) // ' = ', CH % INFO_1
  print *, trim ( CH % LABEL ( 4 ) ) // ' = ', CH % INFO_2
  print *, trim ( CH % LABEL ( 5 ) ) // ' = ', CH % INFO_3
  print *, trim ( CH % LABEL ( 6 ) ) // ' = ', CH % INFO_4
  print *, trim ( CH % LABEL ( 7 ) ) // ' = ', CH % INFO_5
  print *, trim ( CH % LABEL ( 8 ) ) // ' = ', CH % INFO_6
  print *, trim ( CH % LABEL ( 9 ) ) // ' = ', CH % INFO_7

end program ConsoleHeader_Form_Test
