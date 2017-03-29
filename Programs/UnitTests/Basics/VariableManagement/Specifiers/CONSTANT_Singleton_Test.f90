program CONSTANT_Singleton_Test

  use CONSTANT_Singleton
  
  implicit none

  print*
  print*, 'PI                = ', CONSTANT % PI
  print*, 'SPEED_OF_LIGHT    = ', CONSTANT % SPEED_OF_LIGHT
  print*, 'PERMEABILITY      = ', CONSTANT % PERMEABILITY
  print*, 'GRAVITATIONAL     = ', CONSTANT % GRAVITATIONAL
  print*, 'BOLTZMANN         = ', CONSTANT % BOLTZMANN
  print*, 'METER             = ', CONSTANT % METER
  print*, 'SECOND            = ', CONSTANT % SECOND
  print*, 'KILOGRAM          = ', CONSTANT % KILOGRAM
  print*, 'KELVIN            = ', CONSTANT % KELVIN
  print*, 'AMPERE            = ', CONSTANT % AMPERE
  print*, 'ELECTRON_VOLT     = ', CONSTANT % ELECTRON_VOLT
  print*, 'PLANCK_REDUCED    = ', CONSTANT % PLANCK_REDUCED
  print*, 'RADIATION         = ', CONSTANT % RADIATION
  print*, 'ATOMIC_MASS_UNIT  = ', CONSTANT % ATOMIC_MASS_UNIT
  print*, 'ASTRONOMICAL_UNIT = ', CONSTANT % ASTRONOMICAL_UNIT
  print*, 'SOLAR_MASS        = ', CONSTANT % SOLAR_MASS

end program CONSTANT_Singleton_Test
