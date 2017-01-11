!-- CONSTANT_Singleton defines physical and astrophysical constants

module CONSTANT_Singleton

  use KIND_DEFAULT_Singleton
  
  implicit none
  private
  
    real ( KDR ), private, parameter :: &
      !-- Mathematical
      PI = acos ( -1.0_KDR ), &
      !-- Physical 
      !   http://pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf
      SPEED_OF_LIGHT_MKS = 2.99792458e+8_KDR, &
      PLANCK_REDUCED_MKS = 1.054571800e-34_KDR, &
      PERMEABILITY_MKS   = 4.0e-7_KDR * PI, & 
      GRAVITATIONAL_MKS  = 6.67408e-11_KDR, &
      AVOGADRO_MKS       = 6.022140857e+23_KDR, &
      BOLTZMANN_MKS      = 1.38064852e-23_KDR, &
      ELECTRON_VOLT_MKS  = 1.6021766208e-19_KDR, &
      !-- Astrophysical 
      !   http://pdg.lbl.gov/2016/reviews
      !          /rpp2016-rev-astrophysical-constants.pdf
      ASTRONOMICAL_UNIT_MKS = 1.49597870700e+11_KDR, &
      SOLAR_MASS_MKS        = 1.98848e+30_KDR

  type, public :: ConstantSingleton
    real ( KDR ) :: &
      !-- Mathematical
      PI = PI, &
      !-- Physical 
      SPEED_OF_LIGHT_MKS = SPEED_OF_LIGHT_MKS, &
      PLANCK_REDUCED_MKS = PLANCK_REDUCED_MKS, &
      PERMEABILITY_MKS   = PERMEABILITY_MKS, & 
      GRAVITATIONAL_MKS  = GRAVITATIONAL_MKS, &
      AVOGADRO_MKS       = AVOGADRO_MKS, &
      BOLTZMANN_MKS      = BOLTZMANN_MKS, &
      ELECTRON_VOLT_MKS  = ELECTRON_VOLT_MKS, &
      !-- Astrophysical
      ASTRONOMICAL_UNIT_MKS = ASTRONOMICAL_UNIT_MKS, &
      SOLAR_MASS_MKS        = SOLAR_MASS_MKS, &
      !-- GenASiS
      SPEED_OF_LIGHT = 1.0_KDR, &
      PERMEABILITY   = 1.0_KDR, &
      GRAVITATIONAL  = 1.0_KDR, &
      BOLTZMANN      = 1.0_KDR, &
      PLANCK_REDUCED = PLANCK_REDUCED_MKS * GRAVITATIONAL_MKS &
                       / SPEED_OF_LIGHT_MKS ** 3
  end type ConstantSingleton
  
  type ( ConstantSingleton ), public, parameter :: &
    CONSTANT = ConstantSingleton ( )

end module CONSTANT_Singleton
