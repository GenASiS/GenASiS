!-- UNIT_Singleton instantiates MeasuredValueForm objects to a commonly used 
!   units of measure in GenASiS internal unit (meter).

module UNIT_Singleton

  use KIND_DEFAULT_Singleton
  use MeasuredValue_Form
  use CONSTANT_Singleton, &
        C => CONSTANT
  
  implicit none
  private
  
  type, public :: UnitSingleton
    type ( MeasuredValueForm ) :: &  !-- Identity
      IDENTITY
    type ( MeasuredValueForm ) :: &  !-- Length
      METER, &
      CENTIMETER, &
      FEMTOMETER, &
      KILOMETER, &
      PARSEC, &
      GIGAPARSEC, &
      ANGSTROM
    type ( MeasuredValueForm ) :: &  !-- Angle
      RADIAN
    type ( MeasuredValueForm ) :: &  !-- Time
      SECOND, &
      MILLISECOND, &
      FEMTOSECOND
    type ( MeasuredValueForm ) :: &  !-- Frequency
      HERTZ, &
      KILOHERTZ
    type ( MeasuredValueForm ) :: &  !-- Mass
      KILOGRAM, &
      GRAM, &
      ATOMIC_MASS_UNIT, &
      SOLAR_MASS
    type ( MeasuredValueForm ) :: &  !-- Speed
      SPEED_OF_LIGHT
    type ( MeasuredValueForm ) :: &  !-- Energy
      JOULE, &
      ERG, &
      BETHE, &
      ELECTRON_VOLT, &
      MEV
    type ( MeasuredValueForm ) :: &  !-- Force
      NEWTON, &
      DYNE
    type ( MeasuredValueForm ) :: &  !-- Pressure
      PASCAL, &
      BARYE
    type ( MeasuredValueForm ) :: &  !-- Temperature
      KELVIN
    type ( MeasuredValueForm ) :: &  !-- Magnetic current
      AMPERE
    type ( MeasuredValueForm ) :: &  !-- Magnetic field
      TESLA, &
      GAUSS
    type ( MeasuredValueForm ) :: &  !-- Number density
      NUMBER_DENSITY_ANGSTROM
    type ( MeasuredValueForm ) :: &  !-- Mass density
      MASS_DENSITY_CGS
    type ( MeasuredValueForm ) :: &  !-- Energy/length conversion
      HBAR_C
    type ( MeasuredValueForm ) :: &  !-- Computer resources
      KILOBYTE, &
      WALL_TIME
  contains
    procedure, public, nopass :: &
      Initialize
    procedure, public, nopass :: &
      Select
  end type UnitSingleton 
  
  type ( UnitSingleton ), public, protected, save, target :: &
    UNIT
  
contains

  
  subroutine Initialize ( )
    
    associate ( U => UNIT )
    
    !-- Identity
    call U % IDENTITY % Initialize ( '', '', 1.0_KDR )
    
    !-- Length
    call U % METER % Initialize ( 'm', 'm', 1.0_KDR )
    call U % CENTIMETER % Initialize ( 1.0e-2_KDR * U % METER, 'cm' )
    call U % FEMTOMETER % Initialize ( 1.0e-15_KDR * U % METER, 'fm' )
    call U % KILOMETER % Initialize  ( 1.0e+3_KDR * U % METER, 'km' )
    call U % PARSEC % Initialize &
           ( ( C % ASTRONOMICAL_UNIT_MKS * U % METER ) &
             / tan ( 2.0_KDR * C % PI / ( 360.0_KDR * 60.0_KDR * 60.0_KDR ) ), &
            'pc' )
    call U % GIGAPARSEC % Initialize ( 1.0e+9_KDR * U % PARSEC, 'Gpc' )
    call U % ANGSTROM % Initialize  ( 1.0e-10_KDR * U % METER, 'A' )
    
    !-- Angle
    call U % RADIAN % Initialize ( UNIT % IDENTITY, 'rad' )

    !-- Time
    call U % SECOND % Initialize &
           ( C % SPEED_OF_LIGHT_MKS / C % SPEED_OF_LIGHT * U % METER, 's' )
    call U % MILLISECOND % Initialize ( 1.0e-3_KDR * U % SECOND, 'ms' )
    call U % FEMTOSECOND % Initialize ( 1.0e-15_KDR * U % SECOND, 'fs' )

    !-- Frequency
    call U % HERTZ % Initialize ( 1.0_KDR / U % SECOND, 'Hz' )
    call U % KILOHERTZ % Initialize ( 1.0e+3_KDR * U % HERTZ, 'kHz' )
    
    !-- Mass
    call U % KILOGRAM % Initialize &
           ( C % GRAVITATIONAL_MKS / C % GRAVITATIONAL &
             * U % METER ** 3 / U % SECOND ** 2, 'kg' )
    call U % GRAM % Initialize ( 1.0e-3_KDR * U % KILOGRAM, 'g' )
    call U % ATOMIC_MASS_UNIT % Initialize &
           ( 1.0_KDR / C % AVOGADRO_MKS * U % GRAM, 'amu' )
    call U % SOLAR_MASS % Initialize &
           ( C % SOLAR_MASS_MKS * U % KILOGRAM, 'M_sun' )
           
    !-- Speed
    call U % SPEED_OF_LIGHT % Initialize ( 'c', '', C % SPEED_OF_LIGHT )

    !-- Energy
    call U % JOULE % Initialize &
           ( U % KILOGRAM * ( U % METER / U % SECOND ) ** 2, 'J' )
    call U % ERG % Initialize &
           ( U % GRAM * ( U % CENTIMETER / U % SECOND ) ** 2, 'erg' )
    call U % BETHE % Initialize ( 1.0e51_KDR * U % ERG, 'Bethe' )
    call U % ELECTRON_VOLT % Initialize &
           ( C % ELECTRON_VOLT_MKS * U % JOULE, 'eV' )
    call U % MEV % Initialize ( 1.0e6_KDR * U % ELECTRON_VOLT, 'MeV' )

    !-- Force
    call U % NEWTON % Initialize ( U % JOULE / U % METER, 'N' )
    call U % DYNE   % Initialize ( U % ERG / U % CENTIMETER, 'dyn' )

    !-- Pressure
    call U % PASCAL % Initialize ( U % NEWTON / U % METER ** 2, 'Pa' )
    call U % BARYE  % Initialize ( U % DYNE / U % CENTIMETER ** 2, 'Ba' )

    !-- Temperature
    call U % KELVIN % Initialize &
           ( C % BOLTZMANN_MKS / C % BOLTZMANN * U % KILOGRAM &
             * U % METER**2 / U % SECOND ** 2, 'K' )

    !-- Magnetic current
    call U % AMPERE % Initialize &
           ( ( ( C % PERMEABILITY_MKS / C % PERMEABILITY ) &
               * U % NEWTON ) ** 0.5_KDR , 'A' )

    !-- Magnetic field
    call U % TESLA % Initialize &
           ( U % NEWTON / ( U % AMPERE * U % METER ), 'T' )
    call U % GAUSS % Initialize ( 1.0e-4_KDR * U % TESLA, 'G' )

    !-- Number density
    call U % NUMBER_DENSITY_ANGSTROM % Initialize &
           ( 1 / U % ANGSTROM ** 3, 'A^-3' ) 

    !-- Mass density
    call U % MASS_DENSITY_CGS % Initialize &
           ( U % GRAM / U % CENTIMETER ** 3, 'g cm^-3' ) 

    !-- Energy/length conversion
    call U % HBAR_C % Initialize &
           ( '(hBar_c)', 'm^2', C % PLANCK_REDUCED * C % SPEED_OF_LIGHT )  

    !-- Computer resources
    call U % KILOBYTE % Initialize ( 'kB', 'kB', 1.0_KDR )
    call U % WALL_TIME % Initialize ( 's', 's', 1.0_KDR )

    end associate

  end subroutine Initialize
  
  
  subroutine Select ( Selector, Result )
    
    character ( * ), intent ( in ) :: &
      Selector
    type ( MeasuredValueForm ), intent ( out ) :: &
      Result

    select case ( trim ( Selector ) )
    case ( 'METER' )
      Result = UNIT % METER
    case ( 'CENTIMETER' )
      Result = UNIT % CENTIMETER
    case ( 'FEMTOMETER' )    
      Result = UNIT % FEMTOMETER
    case ( 'KILOMETER' )
      Result = UNIT % KILOMETER
    case ( 'PARSEC' )   
      Result = UNIT % PARSEC
    case ( 'GIGAPARSEC' )  
      Result = UNIT % GIGAPARSEC
    case ( 'ANGSTROM' )  
      Result = UNIT % ANGSTROM
    case ( 'RADIAN' )
      Result = UNIT % RADIAN
    case ( 'SECOND' )
      Result = UNIT % SECOND
    case ( 'MILLISECOND' ) 
      Result = UNIT % MILLISECOND
    case ( 'FEMTOSECOND' ) 
      Result = UNIT % FEMTOSECOND
    case ( 'HERTZ' )
      Result = UNIT % HERTZ
    case ( 'KILOHERTZ' )  
      Result = UNIT % KILOHERTZ
    case ( 'KILOGRAM' )      
      Result = UNIT % KILOGRAM
    case ( 'GRAM' )
      Result = UNIT % GRAM
    case ( 'ATOMIC_MASS_UNIT' )
      Result = UNIT % ATOMIC_MASS_UNIT
    case ( 'SOLAR_MASS' )
      Result = UNIT % SOLAR_MASS
    case ( 'SPEED_OF_LIGHT' ) 
      Result = UNIT % SPEED_OF_LIGHT
    case ( 'JOULE' )
      Result = UNIT % JOULE    
    case ( 'ERG' )
      Result = UNIT % ERG
    case ( 'BETHE' )
      Result = UNIT % BETHE   
    case ( 'ELECTRON_VOLT' )
      Result = UNIT % ELECTRON_VOLT
    case ( 'MEV' )
      Result = UNIT % MEV
    case ( 'NEWTON' )
      Result = UNIT % NEWTON
    case ( 'DYNE' )
      Result = UNIT % DYNE
    case ( 'PASCAL' )
      Result = UNIT % PASCAL
    case ( 'BARYE' )
      Result = UNIT % BARYE
    case ( 'KELVIN' )
      Result = UNIT % KELVIN
    case ( 'AMPERE' )
      Result = UNIT % AMPERE
    case ( 'TESLA' )
      Result = UNIT % TESLA
    case ( 'GAUSS' )
      Result = UNIT % GAUSS
    case ( 'NUMBER_DENSITY_ANGSTROM' )
      Result = UNIT % NUMBER_DENSITY_ANGSTROM
    case ( 'MASS_DENSITY_CGS' )
      Result = UNIT % MASS_DENSITY_CGS
    case ( 'HBAR_C' )
      Result = UNIT % HBAR_C
    case ( 'KILOBYTE' )
      Result = UNIT % KILOBYTE
    case ( 'WALL_TIME' )
      Result = UNIT % WALL_TIME
    case default
      Result % Number = 1.0_KDR
      Result % Unit   = 'Undefined'
    end select

  end subroutine Select
  
 
end module UNIT_Singleton
