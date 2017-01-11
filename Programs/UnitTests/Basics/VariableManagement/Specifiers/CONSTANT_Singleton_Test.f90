program CONSTANT_Singleton_Test

  use CONSTANT_Singleton
  
  implicit none

  print*
  print*, 'Mathematical Constants'
  print*, 'PI = ', CONSTANT % PI

  print*
  print*, 'Physical Constants'  
  print*, 'SPEED_OF_LIGHT_MKS = ', CONSTANT % SPEED_OF_LIGHT_MKS
  print*, 'PLANCK_REDUCED_MKS = ', CONSTANT % PLANCK_REDUCED_MKS
  print*, 'PERMEABILITY_MKS   = ', CONSTANT % PERMEABILITY_MKS
  print*, 'GRAVITATIONAL_MKS  = ', CONSTANT % GRAVITATIONAL_MKS
  print*, 'AVOGADRO_MKS       = ', CONSTANT % AVOGADRO_MKS
  print*, 'BOLTZMANN_MKS      = ', CONSTANT % BOLTZMANN_MKS
  print*, 'ELECTRON_VOLT_MKS  = ', CONSTANT % ELECTRON_VOLT_MKS

  print*
  print*, 'Astrophysical Constants'
  print*, 'ASTRONOMICAL_UNIT_MKS = ', CONSTANT % ASTRONOMICAL_UNIT_MKS
  print*, 'SOLAR_MASS_MKS        = ', CONSTANT % SOLAR_MASS_MKS

  print*
  print*, 'GenASiS Constants'
  print*, 'SPEED_OF_LIGHT = ', CONSTANT % SPEED_OF_LIGHT
  print*, 'PERMEABILITY   = ', CONSTANT % PERMEABILITY
  print*, 'GRAVITATIONAL  = ', CONSTANT % GRAVITATIONAL
  print*, 'BOLTZMANN      = ', CONSTANT % BOLTZMANN
  print*, 'PLANCK_REDUCED = ', CONSTANT % PLANCK_REDUCED

end program CONSTANT_Singleton_Test
