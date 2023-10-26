module Units_F__Form
  
  !-- Units_Fluid__Form

  use Basics

  implicit none
  private

  type, public :: Units_F_Form
    !-- Phase space 
    type ( QuantityForm ) :: &
      Time, &
      Length, &
      SqrtDet_M  !-- SquareRoot_Determinant_Metric
!-- FIXME: GCC 11.3 doesn't like hardwired dimensionality
!    type ( MeasuredValueForm ), dimension ( 3 ) :: &
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Coordinate_PS, &  !-- Coordinate_PositionSpace
      Coordinate_MS     !-- Coordinate_MomentumSpace
    !-- Local
    type ( QuantityForm ) :: &
      BaryonMass, &
      NumberDensity, &
      MassDensity, &
      EnergyDensity, &
      Temperature
!-- FIXME: GCC 11.3 doesn't like hardwired dimensionality
!    type ( MeasuredValueForm ), dimension ( 3 ) :: &
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Velocity_U, &
      MomentumDensity_D
    !-- Global
    type ( QuantityForm ) :: &
      Number, &
      Mass, &
      Energy, &
      Momentum, &
      AngularMomentum
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_U
    final :: &
      Finalize
  end type Units_F_Form

contains


  subroutine Initialize ( U, TypeOption, CoordinateSystemOption )

    class ( Units_F_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ), optional :: &
      TypeOption, &
      CoordinateSystemOption

    allocate ( U % Coordinate_PS ( 3 ) )
    allocate ( U % Coordinate_MS ( 3 ) )
    allocate ( U % Velocity_U ( 3 ) )
    allocate ( U % MomentumDensity_D ( 3 ) )

    if ( present ( TypeOption ) ) then
      select case ( trim ( TypeOption ) )
      case ( 'MKS' )
        !-- Phase space 
        U % Time                 =  UNIT % SECOND
        U % Length               =  UNIT % METER
        U % SqrtDet_M            =  UNIT % IDENTITY
        U % Coordinate_PS ( 1 )  =  UNIT % METER
        U % Coordinate_PS ( 2 )  =  UNIT % METER
        U % Coordinate_PS ( 3 )  =  UNIT % METER
        U % Coordinate_MS ( 1 )  =  UNIT % JOULE
        U % Coordinate_MS ( 2 )  =  UNIT % RADIAN
        U % Coordinate_MS ( 3 )  =  UNIT % RADIAN
        !-- Local
        U % BaryonMass               =  UNIT % KILOGRAM
        U % NumberDensity            =  UNIT % NUMBER_DENSITY_MKS
        U % MassDensity              =  UNIT % MASS_DENSITY_MKS
        U % EnergyDensity            =  UNIT % ENERGY_DENSITY_MKS
        U % Temperature              =  UNIT % KELVIN
        U % Velocity_U ( 1 )         =  UNIT % SPEED_MKS
        U % Velocity_U ( 2 )         =  UNIT % SPEED_MKS
        U % Velocity_U ( 3 )         =  UNIT % SPEED_MKS
        U % MomentumDensity_D ( 1 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  UNIT % SPEED_MKS
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  UNIT % SPEED_MKS
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  UNIT % SPEED_MKS
        !-- Global
        U % Number           =  UNIT % MOLE
        U % Mass             =  UNIT % KILOGRAM
        U % Energy           =  UNIT % JOULE
        U % Momentum         =  UNIT % KILOGRAM  *  UNIT % SPEED_MKS
        U % AngularMomentum  =  U % Momentum  *  UNIT % METER
      case ( 'CGS' )
        !-- Phase space 
        U % Time                 =  UNIT % SECOND
        U % Length               =  UNIT % CENTIMETER
        U % SqrtDet_M            =  UNIT % IDENTITY
        U % Coordinate_PS ( 1 )  =  UNIT % CENTIMETER
        U % Coordinate_PS ( 2 )  =  UNIT % CENTIMETER
        U % Coordinate_PS ( 3 )  =  UNIT % CENTIMETER
        U % Coordinate_MS ( 1 )  =  UNIT % ERG
        U % Coordinate_MS ( 2 )  =  UNIT % RADIAN
        U % Coordinate_MS ( 3 )  =  UNIT % RADIAN
        !-- Local
        U % BaryonMass               =  UNIT % GRAM
        U % NumberDensity            =  UNIT % NUMBER_DENSITY_CGS
        U % MassDensity              =  UNIT % MASS_DENSITY_CGS
        U % EnergyDensity            =  UNIT % ENERGY_DENSITY_CGS
        U % Temperature              =  UNIT % KELVIN
        U % Velocity_U ( 1 )         =  UNIT % SPEED_CGS
        U % Velocity_U ( 2 )         =  UNIT % SPEED_CGS
        U % Velocity_U ( 3 )         =  UNIT % SPEED_CGS
        U % MomentumDensity_D ( 1 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  UNIT % SPEED_CGS
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  UNIT % SPEED_CGS
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  UNIT % SPEED_CGS
        !-- Global
        U % Number           =  UNIT % MOLE
        U % Mass             =  UNIT % GRAM
        U % Energy           =  UNIT % ERG
        U % Momentum         =  UNIT % GRAM  *  UNIT % SPEED_CGS
        U % AngularMomentum  =  U % Momentum  *  UNIT % CENTIMETER
      case ( 'ASTROPHYSICS' )
        !-- Phase space 
        U % Time                 =  UNIT % SECOND
        U % Length               =  UNIT % KILOMETER
        U % SqrtDet_M            =  UNIT % IDENTITY
        U % Coordinate_PS ( 1 )  =  UNIT % KILOMETER
        U % Coordinate_PS ( 2 )  =  UNIT % RADIAN
        U % Coordinate_PS ( 3 )  =  UNIT % RADIAN
        U % Coordinate_MS ( 1 )  =  UNIT % MEGA_ELECTRON_VOLT
        U % Coordinate_MS ( 2 )  =  UNIT % RADIAN
        U % Coordinate_MS ( 3 )  =  UNIT % RADIAN
        !-- Local
        U % BaryonMass               =  UNIT % ATOMIC_MASS_UNIT
        U % NumberDensity            =  UNIT % NUMBER_DENSITY_NUCLEAR
        U % MassDensity              =  UNIT % MASS_DENSITY_CGS
        U % EnergyDensity            =  UNIT % ENERGY_DENSITY_NUCLEAR
        U % Temperature              =  UNIT % MEGA_ELECTRON_VOLT
        U % Velocity_U ( 1 )         =  UNIT % KILOMETER  /  UNIT % SECOND
        U % Velocity_U ( 2 )         =  UNIT % RADIAN     /  UNIT % SECOND
        U % Velocity_U ( 3 )         =  UNIT % RADIAN     /  UNIT % SECOND
        U % MomentumDensity_D ( 1 )  =  UNIT % ENERGY_DENSITY_NUCLEAR &
                                        /  UNIT % SPEED_OF_LIGHT
        U % MomentumDensity_D ( 2 )  =  UNIT % ENERGY_DENSITY_NUCLEAR  &
                                        /  UNIT % SPEED_OF_LIGHT
        U % MomentumDensity_D ( 2 )  =  UNIT % ENERGY_DENSITY_NUCLEAR &
                                        /  UNIT % SPEED_OF_LIGHT
        !-- Global
        U % Number           =  UNIT % SOLAR_BARYON_NUMBER
        U % Mass             =  UNIT % SOLAR_MASS
        U % Energy           =  UNIT % ENERGY_SOLAR_MASS
        U % Momentum         =  UNIT % MOMENTUM_SOLAR_MASS
        U % AngularMomentum  =  UNIT % SOLAR_KERR_PARAMETER
      case default
        call Show ( 'Type not recognized', CONSOLE % ERROR )
        call Show ( 'Units_F__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- TypeOption
    end if !-- TypeOption

    if ( present ( CoordinateSystemOption ) ) then
      select case ( trim ( CoordinateSystemOption ) )
      case ( 'RECTANGULAR' )
        !-- Leave defaults
!      case ( 'CYLINDRICAL' )
!      case ( 'SPHERICAL' )
      case default
        call Show ( 'CoordinateSystem not recognized', CONSOLE % ERROR )
        call Show ( 'Units_F__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- CoordinateSystemOption
    end if !-- CoordinateSystemOption 

  end subroutine Initialize


  subroutine Show_U ( U, IgnorabilityOption )

    class ( Units_F_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call Show ( 'Units -- Phase space', IgnorabilityOption )
    call Show ( U % Time,          'Time',          IgnorabilityOption )
    call Show ( U % Length,        'Length',        IgnorabilityOption )
    call Show ( U % SqrtDet_M,     'SqrtDet_M',     IgnorabilityOption )
    call Show ( U % Coordinate_PS, 'Coordinate_PS', IgnorabilityOption )
    call Show ( U % Coordinate_MS, 'Coordinate_MS', IgnorabilityOption )

    call Show ( 'Units -- Local' )
    call Show ( U % BaryonMass,        'BaryonMass',    IgnorabilityOption )
    call Show ( U % NumberDensity,     'NumberDensity', IgnorabilityOption )
    call Show ( U % MassDensity,       'MassDensity',   IgnorabilityOption )
    call Show ( U % EnergyDensity,     'EnergyDensity', IgnorabilityOption )
    call Show ( U % Temperature,       'Temperature',   IgnorabilityOption )
    call Show ( U % Velocity_U,        'Velocity_U',    IgnorabilityOption )
    call Show ( U % MomentumDensity_D, 'MomentumDensity_D', &
                                                        IgnorabilityOption )

    call Show ( 'Units -- Global' )
    call Show ( U % Number,          'Number', IgnorabilityOption )
    call Show ( U % Mass,            'Mass', IgnorabilityOption )
    call Show ( U % Energy,          'Energy', IgnorabilityOption )
    call Show ( U % Momentum,        'Momentum', IgnorabilityOption )
    call Show ( U % AngularMomentum, 'AngularMomentum', IgnorabilityOption )

  end subroutine Show_U


  impure elemental subroutine Finalize ( U )

    type ( Units_F_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % MomentumDensity_D ) ) &
      deallocate ( U % MomentumDensity_D )
    if ( allocated ( U % Velocity_U ) ) &
      deallocate ( U % Velocity_U )
    if ( allocated ( U % Coordinate_MS ) ) &
      deallocate ( U % Coordinate_MS )
    if ( allocated ( U % Coordinate_PS ) ) &
      deallocate ( U % Coordinate_PS )

  end subroutine Finalize


end module Units_F__Form
