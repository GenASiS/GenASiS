!-- UNIT_Singleton instantiates MeasuredValueForm objects to a commonly used 
!   units of measure in GenASiS internal unit (meter).

module UNIT_Singleton

  use KIND_DEFAULT_Singleton
  use KIND_BIG_Singleton, &
        KB => KIND_BIG
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
      ASTRONOMICAL_UNIT, &
      PARSEC, &
      GIGAPARSEC, &
      ANGSTROM
    type ( MeasuredValueForm ) :: &  !-- Mass
      KILOGRAM, &
      GRAM, &
      ATOMIC_MASS_UNIT, &
      SOLAR_MASS
    type ( MeasuredValueForm ) :: &  !-- Time
      SECOND, &
      MILLISECOND, &
      FEMTOSECOND
    type ( MeasuredValueForm ) :: &  !-- Magnetic current
      AMPERE
    type ( MeasuredValueForm ) :: &  !-- Temperature
      KELVIN
    type ( MeasuredValueForm ) :: &  !-- Amount of substance
      MOLE, &
      SOLAR_BARYON_NUMBER
    type ( MeasuredValueForm ) :: &  !-- Angle
      RADIAN
    type ( MeasuredValueForm ) :: &  !-- Frequency
      HERTZ, &
      KILOHERTZ
    type ( MeasuredValueForm ) :: &  !-- Speed
      SPEED_OF_LIGHT
    type ( MeasuredValueForm ) :: &  !-- Force
      NEWTON, &
      DYNE
    type ( MeasuredValueForm ) :: &  !-- Pressure
      PASCAL, &
      BARYE
    type ( MeasuredValueForm ) :: &  !-- Magnetic field
      TESLA, &
      GAUSS
    type ( MeasuredValueForm ) :: &  !-- Electric potential
      VOLT
    type ( MeasuredValueForm ) :: &  !-- Energy
      JOULE, &
      ERG, &
      ELECTRON_VOLT, &
      MEGA_ELECTRON_VOLT, &
      BETHE
    type ( MeasuredValueForm ) :: &  !-- Angular momentum
      SOLAR_KERR_PARAMETER
    type ( MeasuredValueForm ) :: &  !-- Entropy per baryon
      BOLTZMANN
    type ( MeasuredValueForm ) :: &  !-- Energy/length conversion
      HBAR_C
    type ( MeasuredValueForm ) :: &  !-- Number density
      NUMBER_DENSITY_ANGSTROM, &
      NUMBER_DENSITY_MEV_HBAR_C
    type ( MeasuredValueForm ) :: &  !-- Mass density
      MASS_DENSITY_CGS
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

    character ( 6 ) :: &
      MeV_Minus_1 = 'MeV^-1'
    character ( 5, KBCH ) :: &
      MeV_Minus_1_KBCH &
        = KBCH_'MeV' // char ( KB % SUPERSCRIPT_MINUS, KBCH ) &
                     // char ( KB % SUPERSCRIPT_1, KBCH )

contains

  
  subroutine Initialize ( )
    
    associate ( U => UNIT )

    !-- Identity
    call U % IDENTITY % Initialize &
           ( '', '', 1.0_KDR )
    
    !-- Length
    if ( KBCH > KDCH ) then
      call U % METER % Initialize_UCS &
             ( KBCH_'m', MeV_Minus_1_KBCH, 'm', MeV_Minus_1, C % METER )
    else
      call U % METER % Initialize &
             ( 'm', MeV_Minus_1, C % METER )
    end if
    call U % CENTIMETER % Initialize &
           ( 1.0e-2_KDR * U % METER, 'cm' )
    call U % FEMTOMETER % Initialize &
           ( 1.0e-15_KDR * U % METER, 'fm' )
    call U % KILOMETER % Initialize &
           ( 1.0e+3_KDR * U % METER, 'km' )
    call U % ASTRONOMICAL_UNIT % Initialize &
           ( 'au', MeV_Minus_1, C % ASTRONOMICAL_UNIT )
    call U % PARSEC % Initialize &
           ( 'pc', MeV_Minus_1, C % PARSEC )
    call U % GIGAPARSEC % Initialize &
           ( 1.0e+9_KDR * U % PARSEC, 'Gpc' )
    if ( KBCH > KDCH ) then
      call U % ANGSTROM % Initialize_UCS &
             ( 1.0e-10_KDR * U % METER, &
               char ( KB % CAPITAL_A_RING, KBCH ), 'Ang' )
    else
      call U % ANGSTROM % Initialize &
             ( 1.0e-10_KDR * U % METER, 'Ang' )
    end if

    !-- Mass
    call U % KILOGRAM % Initialize &
           ( 'kg', 'MeV', C % KILOGRAM )
    call U % GRAM % Initialize &
           ( 1.0e-3_KDR * U % KILOGRAM, 'g' )
    call U % ATOMIC_MASS_UNIT % Initialize &
           ( 'u', 'MeV', C % ATOMIC_MASS_UNIT )
    if ( KBCH > KDCH ) then
      call U % SOLAR_MASS % Initialize_UCS &
             ( KBCH_'M' // char ( KB % CIRCLE_DOT, KBCH ), KBCH_'MeV', &
               'M_Sun', 'MeV', C % SOLAR_MASS )
    else
      call U % SOLAR_MASS % Initialize &
             ( 'M_Sun', 'MeV', C % SOLAR_MASS )
    end if
           
    !-- Time
    if ( KBCH > KDCH ) then
      call U % SECOND % Initialize_UCS &
             ( KBCH_'s', MeV_Minus_1_KBCH, 's', MeV_Minus_1, C % SECOND )
    else
      call U % SECOND % Initialize &
             ( 's', MeV_Minus_1, C % SECOND )
    end if
    call U % MILLISECOND % Initialize &
           ( 1.0e-3_KDR * U % SECOND, 'ms' )
    call U % FEMTOSECOND % Initialize &
           ( 1.0e-15_KDR * U % SECOND, 'fs' )

    !-- Magnetic current
    call U % AMPERE % Initialize &
           ( 'A', 'MeV', C % AMPERE )

    !-- Temperature
    call U % KELVIN % Initialize &
           ( 'K', 'MeV', C % KELVIN )

    !-- Amount of substance
    call U % MOLE % Initialize &
           ( 'mol', '', C % MOLE )
    if ( KBCH > KDCH ) then
      call U % SOLAR_BARYON_NUMBER % Initialize_UCS &
             ( KBCH_'N' // char ( KB % CIRCLE_DOT, KBCH ), KBCH_'', &
               'N_Sun', '', C % SOLAR_BARYON_NUMBER )
    else
      call U % SOLAR_BARYON_NUMBER % Initialize &
             ( 'N_Sun', '', C % SOLAR_BARYON_NUMBER )
    end if

    !-- Angle
    call U % RADIAN % Initialize &
           ( UNIT % IDENTITY, 'rad' )

    !-- Frequency
    call U % HERTZ % Initialize &
           ( 1.0_KDR / U % SECOND, 'Hz' )
    call U % KILOHERTZ % Initialize &
           ( 1.0e+3_KDR * U % HERTZ, 'kHz' )
    
    !-- Speed
    call U % SPEED_OF_LIGHT % Initialize &
           ( 'c', '', C % SPEED_OF_LIGHT )

    !-- Force
    call U % NEWTON % Initialize &
           ( U % KILOGRAM  *  U % METER  *  U % SECOND ** (-2), 'N' )
    call U % DYNE % Initialize &
           ( U % GRAM * U % CENTIMETER  *  U % SECOND ** (-2), 'dyn' )

    !-- Pressure
    call U % PASCAL % Initialize &
           ( U % NEWTON / U % METER ** 2, 'Pa' )
    call U % BARYE  % Initialize &
           ( U % DYNE / U % CENTIMETER ** 2, 'Ba' )

    !-- Magnetic field
    call U % TESLA % Initialize &
           ( U % NEWTON / ( U % AMPERE * U % METER ), 'T' )
    call U % GAUSS % Initialize &
           ( 1.0e-4_KDR * U % TESLA, 'G' )

    !-- Electric potential
    call U % VOLT % Initialize &
           ( U % NEWTON * U % METER / ( U % AMPERE * U % SECOND ), 'V' )

    !-- Energy
    call U % JOULE % Initialize &
           ( U % NEWTON * U % METER, 'J' )
    call U % ERG % Initialize &
           ( U % DYNE * U % CENTIMETER, 'erg' )
    call U % ELECTRON_VOLT % Initialize &
           ( C % ELECTRON_CHARGE * U % VOLT, 'eV' )
    call U % MEGA_ELECTRON_VOLT % Initialize &
           ( 1.0e6_KDR * U % ELECTRON_VOLT, 'MeV' )
    call U % BETHE % Initialize &
           ( 1.0e51_KDR * U % ERG, 'B' )

    !-- Angular momentum
    if ( KBCH > KDCH ) then
      call U % SOLAR_KERR_PARAMETER % Initialize_UCS &
             ( KBCH_'G M' // char ( KB % CIRCLE_DOT, KBCH ) &
               // char ( KB % SUPERSCRIPT_2, KBCH ) // KBCH_' c' &
               // char ( KB % SUPERSCRIPT_MINUS, KBCH ) &
               // char ( KB % SUPERSCRIPT_1, KBCH ), &
               KBCH_'', &
               'G M_Sun^2 c^-1', '', &
               C % SOLAR_KERR_PARAMETER )
    else
      call U % SOLAR_KERR_PARAMETER % Initialize &
             ( 'G M_Sun^2 c^-1', '', C % SOLAR_KERR_PARAMETER )
    end if

    !-- Entropy per baryon
    call U % BOLTZMANN % Initialize &
           ( 'k', '', C % BOLTZMANN )

    !-- Energy/length conversion
    if ( KBCH > KDCH ) then
      call U % HBAR_C % Initialize_UCS &
             ( KBCH_'(' // char ( KB % SMALL_H_STROKE, KBCH ) &
                 // KBCH_'c)', &
               KBCH_'', '(hBar_c)', '', &
               C % PLANCK_REDUCED * C % SPEED_OF_LIGHT )
    else
      call U % HBAR_C % Initialize &
             ( '(hBar_c)', '', C % PLANCK_REDUCED * C % SPEED_OF_LIGHT )
    end if

    !-- Number density
    U % NUMBER_DENSITY_ANGSTROM &
      =  1 / U % ANGSTROM ** 3
    U % NUMBER_DENSITY_MEV_HBAR_C &
      =  U % MEGA_ELECTRON_VOLT ** 3  /  U % HBAR_C ** 3

    !-- Mass density
    U % MASS_DENSITY_CGS &
      =  U % GRAM  /  U % CENTIMETER ** 3

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
    case ( 'IDENTITY' )
      Result = UNIT % IDENTITY
    case ( 'METER' )
      Result = UNIT % METER
    case ( 'CENTIMETER' )
      Result = UNIT % CENTIMETER
    case ( 'FEMTOMETER' )    
      Result = UNIT % FEMTOMETER
    case ( 'KILOMETER' )
      Result = UNIT % KILOMETER
    case ( 'ASTRONOMICAL_UNIT' )
      Result = UNIT % ASTRONOMICAL_UNIT
    case ( 'PARSEC' )   
      Result = UNIT % PARSEC
    case ( 'GIGAPARSEC' )  
      Result = UNIT % GIGAPARSEC
    case ( 'ANGSTROM' )  
      Result = UNIT % ANGSTROM
    case ( 'KILOGRAM' )      
      Result = UNIT % KILOGRAM
    case ( 'GRAM' )
      Result = UNIT % GRAM
    case ( 'ATOMIC_MASS_UNIT' )
      Result = UNIT % ATOMIC_MASS_UNIT
    case ( 'SOLAR_MASS' )
      Result = UNIT % SOLAR_MASS
    case ( 'SECOND' )
      Result = UNIT % SECOND
    case ( 'MILLISECOND' ) 
      Result = UNIT % MILLISECOND
    case ( 'FEMTOSECOND' ) 
      Result = UNIT % FEMTOSECOND
    case ( 'AMPERE' )
      Result = UNIT % AMPERE
    case ( 'KELVIN' )
      Result = UNIT % KELVIN
    case ( 'MOLE' )
      Result = UNIT % MOLE
    case ( 'SOLAR_BARYON_NUMBER' )
      Result = UNIT % SOLAR_BARYON_NUMBER
    case ( 'RADIAN' )
      Result = UNIT % RADIAN
    case ( 'HERTZ' )
      Result = UNIT % HERTZ
    case ( 'KILOHERTZ' )  
      Result = UNIT % KILOHERTZ
    case ( 'SPEED_OF_LIGHT' ) 
      Result = UNIT % SPEED_OF_LIGHT
    case ( 'NEWTON' )
      Result = UNIT % NEWTON
    case ( 'DYNE' )
      Result = UNIT % DYNE
    case ( 'PASCAL' )
      Result = UNIT % PASCAL
    case ( 'BARYE' )
      Result = UNIT % BARYE
    case ( 'TESLA' )
      Result = UNIT % TESLA
    case ( 'GAUSS' )
      Result = UNIT % GAUSS
    case ( 'VOLT' )
      Result = UNIT % VOLT
    case ( 'JOULE' )
      Result = UNIT % JOULE    
    case ( 'ERG' )
      Result = UNIT % ERG
    case ( 'ELECTRON_VOLT' )
      Result = UNIT % ELECTRON_VOLT
    case ( 'MEGA_ELECTRON_VOLT' )
      Result = UNIT % MEGA_ELECTRON_VOLT
    case ( 'BETHE' )
      Result = UNIT % BETHE   
    case ( 'BOLTZMANN' )
      Result = UNIT % BOLTZMANN
    case ( 'HBAR_C' )
      Result = UNIT % HBAR_C
    case ( 'NUMBER_DENSITY_ANGSTROM' )
      Result = UNIT % NUMBER_DENSITY_ANGSTROM
    case ( 'NUMBER_DENSITY_MEV_HBAR_C' )
      Result = UNIT % NUMBER_DENSITY_MEV_HBAR_C
    case ( 'MASS_DENSITY_CGS' )
      Result = UNIT % MASS_DENSITY_CGS
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
