!-- CONSTANT_Singleton defines physical and astrophysical constants

module CONSTANT_Singleton

  use KIND_DEFAULT_Singleton
  
  implicit none
  private
  
    real ( KDR ), private, parameter :: &
      !-- Mathematical
      PI  =  acos ( -1.0_KDR ), &
      !-- Physical MKS
      !   http://pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf
      SPEED_OF_LIGHT_MKS  =  2.99792458e+8_KDR, &
      PLANCK_REDUCED_MKS  =  1.054571800e-34_KDR, &
      PERMEABILITY_MKS    =  4.0e-7_KDR * PI, & 
      GRAVITATIONAL_MKS   =  6.67408e-11_KDR, &
      AVOGADRO_MKS        =  6.022140857e+23_KDR, &
      BOLTZMANN_MKS       =  1.38064852e-23_KDR, &
      ELECTRON_VOLT_MKS   =  1.6021766208e-19_KDR, &
      !-- Astrophysical MKS
      !   http://pdg.lbl.gov/2016/reviews
      !          /rpp2016-rev-astrophysical-constants.pdf
      ASTRONOMICAL_UNIT_MKS = 1.49597870700e+11_KDR, &
      SOLAR_MASS_MKS        = 1.98848e+30_KDR, &
      !-- GenASiS Fundamental
      SPEED_OF_LIGHT  =  1.0_KDR, &
      PERMEABILITY    =  1.0_KDR, &
      GRAVITATIONAL   =  1.0_KDR, &
      BOLTZMANN       =  1.0_KDR, &
      !-- GenASiS Metric
      METER         =  1.0_KDR, &
      SECOND        =  SPEED_OF_LIGHT_MKS / SPEED_OF_LIGHT &
                       *  METER, &
      KILOGRAM      =  GRAVITATIONAL_MKS / GRAVITATIONAL &
                       *  METER ** 3  /  SECOND ** 2, &
      KELVIN        =  BOLTZMANN_MKS / BOLTZMANN &
                       *  KILOGRAM  *  METER ** 2  /  SECOND ** 2, &
      AMPERE        =  sqrt ( ( PERMEABILITY_MKS / PERMEABILITY ) &
                              *  KILOGRAM  *  METER  /  SECOND ** 2 ), &
      ELECTRON_VOLT =  ELECTRON_VOLT_MKS &
                       *  KILOGRAM  *  METER ** 2  /  SECOND ** 2, &
      !-- GenASiS Derived
      PLANCK_REDUCED    =  PLANCK_REDUCED_MKS  &
                           *  KILOGRAM  *  METER ** 2  /  SECOND, &    
      RADIATION         =  PI ** 2  /  15.0_KDR &
                           *  BOLTZMANN ** 4 &
                           /  ( PLANCK_REDUCED * SPEED_OF_LIGHT ) ** 3, &
      ATOMIC_MASS_UNIT  =  AVOGADRO_MKS ** (-1) &
                           *  1.0e-3_KDR  *  KILOGRAM, &
      ASTRONOMICAL_UNIT =  ASTRONOMICAL_UNIT_MKS &
                           *  METER, &
      SOLAR_MASS        =  SOLAR_MASS_MKS * KILOGRAM

  type, public :: ConstantSingleton
    real ( KDR ) :: &
      PI                 =  PI, &
      SPEED_OF_LIGHT     =  SPEED_OF_LIGHT, &
      PERMEABILITY       =  PERMEABILITY, &
      GRAVITATIONAL      =  GRAVITATIONAL, &
      BOLTZMANN          =  BOLTZMANN, &
      METER              =  METER, &
      SECOND             =  SECOND, &
      KILOGRAM           =  KILOGRAM, &
      KELVIN             =  KELVIN, &
      AMPERE             =  AMPERE, &
      ELECTRON_VOLT      =  ELECTRON_VOLT, &
      PLANCK_REDUCED     =  PLANCK_REDUCED, &
      RADIATION          =  RADIATION, &
      ATOMIC_MASS_UNIT   =  ATOMIC_MASS_UNIT, &
      ASTRONOMICAL_UNIT  =  ASTRONOMICAL_UNIT, &
      SOLAR_MASS         =  SOLAR_MASS
  end type ConstantSingleton
  
  type ( ConstantSingleton ), public, parameter :: &
    CONSTANT = ConstantSingleton ( )

end module CONSTANT_Singleton
