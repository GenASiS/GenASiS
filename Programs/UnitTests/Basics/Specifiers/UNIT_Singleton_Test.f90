program UNIT_Singleton_Test

  use ISO_FORTRAN_ENV
  use KIND_BIG_Singleton
  use Quantity_Form
  use UNIT_Singleton
  
  implicit none
  
  character ( 5 ) :: &
    Encoding

!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  call UNIT % Initialize ( )

  print*
  print*, 'Identity Unit'
  print*, 'IDENTITY            ', UNIT % IDENTITY 

  print*
  print*, 'Length Units'
  print*, 'METER               ', UNIT % METER      
  print*, 'CENTIMETER          ', UNIT % CENTIMETER 
  print*, 'FEMTOMETER          ', UNIT % FEMTOMETER 
  print*, 'KILOMETER           ', UNIT % KILOMETER
  print*, 'ASTRONOMICAL_UNIT   ', UNIT % ASTRONOMICAL_UNIT  
  print*, 'PARSEC              ', UNIT % PARSEC     
  print*, 'GIGAPARSEC          ', UNIT % GIGAPARSEC
  print*, 'ANGSTROM            ', UNIT % ANGSTROM

  print*
  print*, 'Mass Units'  
  print*, 'KILOGRAM            ', UNIT % KILOGRAM   
  print*, 'GRAM                ', UNIT % GRAM       
  print*, 'ATOMIC_MASS_UNIT    ', UNIT % ATOMIC_MASS_UNIT
  print*, 'SOLAR_MASS          ', UNIT % SOLAR_MASS

  print*
  print*, 'Time Units' 
  print*, 'SECOND              ', UNIT % SECOND     
  print*, 'MILLISECOND         ', UNIT % MILLISECOND
  print*, 'FEMTOSECOND         ', UNIT % FEMTOSECOND

  print*
  print*, 'Magnetic Current Units'
  print*, 'AMPERE              ', UNIT % AMPERE

  print*
  print*, 'Temperature Units'
  print*, 'KELVIN              ', UNIT % KELVIN

  print*
  print*, 'Amount of Substance Units'
  print*, 'MOLE                ', UNIT % MOLE
  print*, 'SOLAR_BARYON_NUMBER ', UNIT % SOLAR_BARYON_NUMBER

  print*
  print*, 'Angle Units' 
  print*, 'RADIAN              ', UNIT % RADIAN     

  print*
  print*, 'Frequency Units'
  print*, 'HERTZ               ', UNIT % HERTZ      
  print*, 'KILOHERTZ           ', UNIT % KILOHERTZ

  print*
  print*, 'Speed Units'
  print*, 'SPEED_OF_LIGHT      ', UNIT % SPEED_OF_LIGHT

  print*
  print*, 'Momentum Units'
  print*, 'MOMENTUM_SOLAR_MASS', UNIT % MOMENTUM_SOLAR_MASS

  print*
  print*, 'Force Units'
  print*, 'NEWTON              ', UNIT % NEWTON
  print*, 'DYNE                ', UNIT % DYNE

  print*
  print*, 'Pressure Units'
  print*, 'PASCAL              ', UNIT % PASCAL
  print*, 'BARYE               ', UNIT % BARYE

  print*
  print*, 'Magnetic Field Units'
  print*, 'TESLA               ', UNIT % TESLA
  print*, 'GAUSS               ', UNIT % GAUSS

  print*
  print*, 'Electric Potential Units'
  print*, 'VOLT                ', UNIT % VOLT

  print*
  print*, 'Energy Units'
  print*, 'JOULE               ', UNIT % JOULE
  print*, 'ERG                 ', UNIT % ERG
  print*, 'ELECTRON_VOLT       ', UNIT % ELECTRON_VOLT
  print*, 'MEGA_ELECTRON_VOLT  ', UNIT % MEGA_ELECTRON_VOLT
  print*, 'BETHE               ', UNIT % BETHE
  print*, 'ENERGY_SOLAR_MASS   ', UNIT % ENERGY_SOLAR_MASS

  print*
  print*, 'Angular momentum Units'
  print*, 'SOLAR_KERR_PARAMETER', UNIT % SOLAR_KERR_PARAMETER

  print*
  print*, 'Entropy per baryon Units'
  print*, 'BOLTZMANN           ', UNIT % BOLTZMANN

  print*
  print*, 'Energy/length conversion'
  print*, 'HBAR_C              ', UNIT % HBAR_C

  print*
  print*, 'Number Density Units'
  print*, 'NUMBER_DENSITY_ANGSTROM  ', UNIT % NUMBER_DENSITY_ANGSTROM
  print*, 'NUMBER_DENSITY_NUCLEAR   ', UNIT % NUMBER_DENSITY_NUCLEAR
  print*, 'NUMBER_DENSITY_MEV_HBAR_C', UNIT % NUMBER_DENSITY_MEV_HBAR_C

  print*
  print*, 'Mass Density Units'
  print*, 'MASS_DENSITY_CGS', UNIT % MASS_DENSITY_CGS

  print*
  print*, 'Energy Density Units'
  print*, 'ENERGY_DENSITY_NUCLEAR', UNIT % ENERGY_DENSITY_NUCLEAR

  print*
  print*, 'Computer Resources Units'
  print*, 'KILOBYTE ', UNIT % KILOBYTE
  print*, 'WALL_TIME', UNIT % WALL_TIME

end program UNIT_Singleton_Test
