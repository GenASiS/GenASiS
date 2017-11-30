program CONSTANT_Singleton_Test

  use CONSTANT_Singleton
  
  implicit none

  print*
  print*, 'PI                   = ', CONSTANT % PI
  print*, 'METER                = ', CONSTANT % METER
  print*, 'KILOGRAM             = ', CONSTANT % KILOGRAM
  print*, 'SECOND               = ', CONSTANT % SECOND
  print*, 'AMPERE               = ', CONSTANT % AMPERE
  print*, 'KELVIN               = ', CONSTANT % KELVIN
  print*, 'MOLE                 = ', CONSTANT % MOLE
  print*, 'SPEED_OF_LIGHT       = ', CONSTANT % SPEED_OF_LIGHT
  print*, 'PLANCK_REDUCED       = ', CONSTANT % PLANCK_REDUCED
  print*, 'ELECTRON_CHARGE      = ', CONSTANT % ELECTRON_CHARGE
  print*, 'ELECTRON_MASS        = ', CONSTANT % ELECTRON_MASS
  print*, 'PROTON_MASS          = ', CONSTANT % PROTON_MASS
  print*, 'ATOMIC_MASS_UNIT     = ', CONSTANT % ATOMIC_MASS_UNIT
  print*, 'PERMEABILITY         = ', CONSTANT % PERMEABILITY
  print*, 'GRAVITATIONAL        = ', CONSTANT % GRAVITATIONAL
  print*, 'BOLTZMANN            = ', CONSTANT % BOLTZMANN
  print*, 'STEFAN_BOLTZMANN     = ', CONSTANT % STEFAN_BOLTZMANN
  print*, 'FERMI_COUPLING       = ', CONSTANT % FERMI_COUPLING
  print*, 'ASTRONOMICAL_UNIT    = ', CONSTANT % ASTRONOMICAL_UNIT
  print*, 'PARSEC               = ', CONSTANT % PARSEC
  print*, 'SOLAR_MASS           = ', CONSTANT % SOLAR_MASS
  print*, 'SOLAR_BARYON_NUMBER  = ', CONSTANT % SOLAR_BARYON_NUMBER
  print*, 'SOLAR_KERR_PARAMETER = ', CONSTANT % SOLAR_KERR_PARAMETER
  print*, 'NEUTRON_MASS         = ', CONSTANT % NEUTRON_MASS

end program CONSTANT_Singleton_Test
