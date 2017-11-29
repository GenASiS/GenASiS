!-- CONSTANT_Singleton defines physical and astrophysical constants

module CONSTANT_Singleton

  !-- "Natural units" with Lorentz-Heaviside electron charge, base unit = MeV
  !     ( hBar = c = k = \mu = MeV = 1 )
  !   https://en.wikipedia.org/wiki
  !     /Natural_units#.22Natural_units.22_.28particle_physics_and_cosmology.29

  use KIND_DEFAULT_Singleton
  
  implicit none
  private
  
    real ( KDR ), private, parameter :: &
      !-- Mathematical
      PI  =  acos ( -1.0_KDR ), &
      !-- Physical SI
      !   http://pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf
      SPEED_OF_LIGHT_SI     =  2.99792458e+8_KDR, &
      PLANCK_REDUCED_SI     =  1.054571800e-34_KDR, &
      ELECTRON_CHARGE_SI    =  1.6021766208e-19_KDR, &
      ELECTRON_MASS_SI      =  9.10938356e-31_KDR, &
      PROTON_MASS_SI        =  1.672621898e-27_KDR, &
      ATOMIC_MASS_UNIT_SI   =  1.660539040e-27_KDR, &
      PERMEABILITY_SI       =  4.0e-7_KDR * PI, & 
      GRAVITATIONAL_SI      =  6.67408e-11_KDR, &
      AVOGADRO_SI           =  6.022140857e23_KDR, &
      BOLTZMANN_SI          =  1.38064852e-23_KDR, &
      FERMI_HBAR_C_3_GEV_2  =  1.1663787e-5_KDR, &
      SIN_2_WEINBERG        =  0.23129, &
      !-- Astrophysical SI
      !   http://pdg.lbl.gov/2016/reviews
      !          /rpp2016-rev-astrophysical-constants.pdf
      ASTRONOMICAL_UNIT_SI  =  1.49597870700e+11_KDR, &
      PARSEC_SI             =  3.08567758149e+16_KDR, &
      SOLAR_MASS_SI         =  1.98848e+30_KDR, &
      !-- Particles (that are not in the PDG Physical Constants)  
      !   http://www-pdg.lbl.gov/2016/tables/rpp2016-tab-baryons-N.pdf
      NEUTRON_MASS_AMU        =   1.0086649159_KDR, &
      NEUTRON_AXIAL_COUPLING  =  -1.2723_KDR, &
      !-- GenASiS Fundamental
      SPEED_OF_LIGHT  =  1.0_KDR, &
      PLANCK_REDUCED  =  1.0_KDR, &    
      PERMEABILITY    =  1.0_KDR, &
      BOLTZMANN       =  1.0_KDR, &
      !-- GenASiS Metric
      MEGA_ELECTRON_VOLT  =  1.0_KDR, &
      JOULE               =  1.0e-6_KDR * MEGA_ELECTRON_VOLT &
                             / ELECTRON_CHARGE_SI, & 
      SECOND              =  PLANCK_REDUCED / PLANCK_REDUCED_SI &
                             / JOULE, &
      METER               =  SPEED_OF_LIGHT / SPEED_OF_LIGHT_SI &
                             * SECOND, &
      KILOGRAM            =  JOULE  *  SECOND ** 2  /  METER ** 2, &
      AMPERE              =  sqrt ( ( PERMEABILITY_SI / PERMEABILITY ) &
                                    *  KILOGRAM  *  METER  /  SECOND ** 2 ), &
      KELVIN              =  BOLTZMANN_SI / BOLTZMANN &
                             *  JOULE, &
      MOLE                =  AVOGADRO_SI, &
      !-- GenASiS Derived
      ELECTRON_CHARGE       =  ELECTRON_CHARGE_SI &
                               *  AMPERE * SECOND, &
      ELECTRON_MASS         =  ELECTRON_MASS_SI &
                               *  KILOGRAM, &
      PROTON_MASS           =  PROTON_MASS_SI &
                               *  KILOGRAM, &
      ATOMIC_MASS_UNIT      =  ATOMIC_MASS_UNIT_SI &
                               *  KILOGRAM, &
      GRAVITATIONAL         =  GRAVITATIONAL_SI &
                               *  METER ** 3  /  KILOGRAM  /  SECOND ** 2, &
      STEFAN_BOLTZMANN      =  PI ** 2  /  60.0_KDR &
                               *  BOLTZMANN ** 4 &
                               /  ( PLANCK_REDUCED * SPEED_OF_LIGHT ) ** 3, &
      FERMI_COUPLING        =  FERMI_HBAR_C_3_GEV_2 &
                               *  ( PLANCK_REDUCED * SPEED_OF_LIGHT ) ** 3 &
                               /  ( 1.0e3_KDR * MEGA_ELECTRON_VOLT ) ** 2, &
      ASTRONOMICAL_UNIT     =  ASTRONOMICAL_UNIT_SI &
                               *  METER, &
      PARSEC                =  PARSEC_SI &
                               *  METER, &
      SOLAR_MASS            =  SOLAR_MASS_SI &
                               *  KILOGRAM, &
      SOLAR_BARYON_NUMBER   =  SOLAR_MASS / ATOMIC_MASS_UNIT, &
      SOLAR_KERR_PARAMETER  =  GRAVITATIONAL  *  SOLAR_MASS ** 2  &
                               /  SPEED_OF_LIGHT, & 
      NEUTRON_MASS          =  NEUTRON_MASS_AMU &
                               *  ATOMIC_MASS_UNIT

  type, public :: ConstantSingleton
    real ( KDR ) :: &
      !-- Mathematical
      PI                  =  PI, &
      !-- Metric, SI Base 
      !   https://en.wikipedia.org/wiki/SI_base_unit
      METER     =  METER, &
      KILOGRAM  =  KILOGRAM, &
      SECOND    =  SECOND, &
      AMPERE    =  AMPERE, &
      KELVIN    =  KELVIN, &
      MOLE      =  MOLE, &
      !-- Physical
      SPEED_OF_LIGHT    =  SPEED_OF_LIGHT, &
      PLANCK_REDUCED    =  PLANCK_REDUCED, &
      ELECTRON_CHARGE   =  ELECTRON_CHARGE, &
      ELECTRON_MASS     =  ELECTRON_MASS, &
      PROTON_MASS       =  PROTON_MASS, &
      ATOMIC_MASS_UNIT  =  ATOMIC_MASS_UNIT, &
      PERMEABILITY      =  PERMEABILITY, &
      GRAVITATIONAL     =  GRAVITATIONAL, &
      BOLTZMANN         =  BOLTZMANN, &
      STEFAN_BOLTZMANN  =  STEFAN_BOLTZMANN, &
      FERMI_COUPLING    =  FERMI_COUPLING, &
      SIN_2_WEINBERG    =  SIN_2_WEINBERG, &
      !-- Astrophysical
      ASTRONOMICAL_UNIT       =  ASTRONOMICAL_UNIT, &
      PARSEC                  =  PARSEC, &
      SOLAR_MASS              =  SOLAR_MASS, &
      SOLAR_BARYON_NUMBER     =  SOLAR_BARYON_NUMBER, &
      SOLAR_KERR_PARAMETER    =  SOLAR_KERR_PARAMETER, &
      !-- Particles (that are not in the PDG Physical Constants)
      NEUTRON_MASS            =  NEUTRON_MASS, &
      NEUTRON_AXIAL_COUPLING  =  NEUTRON_AXIAL_COUPLING
  end type ConstantSingleton
  
  type ( ConstantSingleton ), public, parameter :: &
    CONSTANT = ConstantSingleton ( )

end module CONSTANT_Singleton
