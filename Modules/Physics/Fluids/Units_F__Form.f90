module Units_F__Form
  
  !-- Units_Fluid__Form

  use Basics

  implicit none
  private

  type, public :: Units_F_Form
    !-- Spacetime 
    type ( QuantityForm ) :: &
      Time, &
      Length, &
      SqrtDet_M  !-- SquareRoot_Determinant_Metric
!-- FIXME: GCC 11.3 doesn't like hardwired dimensionality
!    type ( MeasuredValueForm ), dimension ( 3 ) :: &
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Coordinate_PS  !-- Coordinate_PositionSpace
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


  subroutine Initialize ( U, CoordinateUnit, TypeOption )

    class ( Units_F_Form ), intent ( inout ) :: &
      U
    type ( QuantityForm ), dimension ( : ), intent ( in ) :: &
      CoordinateUnit
    character ( * ), intent ( in ), optional :: &
      TypeOption

    allocate ( U % Coordinate_PS ( 3 ) )
    allocate ( U % Velocity_U ( 3 ) )
    allocate ( U % MomentumDensity_D ( 3 ) )

    if ( present ( TypeOption ) ) then
      select case ( trim ( TypeOption ) )
      case ( '' )
      case ( 'MKS' )
        !-- Spacetime 
        U % Time           =  UNIT % SECOND
        U % Length         =  UNIT % METER
        U % SqrtDet_M      =  UNIT % IDENTITY
        U % Coordinate_PS  =  CoordinateUnit
        !-- Local
        U % BaryonMass               =  UNIT % KILOGRAM
        U % NumberDensity            =  UNIT % NUMBER_DENSITY_MKS
        U % MassDensity              =  UNIT % MASS_DENSITY_MKS
        U % EnergyDensity            =  UNIT % ENERGY_DENSITY_MKS
        U % Temperature              =  UNIT % KELVIN
        U % Velocity_U ( 1 )         =  U % Coordinate_PS ( 1 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_MKS
        U % Velocity_U ( 2 )         =  U % Coordinate_PS ( 2 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_MKS
        U % Velocity_U ( 3 )         =  U % Coordinate_PS ( 3 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_MKS
        U % MomentumDensity_D ( 1 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 1 ) ** 2  &
                                        *  U % Velocity_U ( 1 )
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 2 ) ** 2  &
                                        *  U % Velocity_U ( 2 )
        U % MomentumDensity_D ( 3 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 3 ) ** 2  &
                                        *  U % Velocity_U ( 3 )
        !-- Global
        U % Number           =  UNIT % MOLE
        U % Mass             =  UNIT % KILOGRAM
        U % Energy           =  UNIT % JOULE
        U % Momentum         =  UNIT % KILOGRAM  *  UNIT % SPEED_MKS
        U % AngularMomentum  =  U % Momentum  *  UNIT % METER
      case ( 'CGS' )
        !-- Spacetime 
        U % Time           =  UNIT % SECOND
        U % Length         =  UNIT % CENTIMETER
        U % SqrtDet_M      =  UNIT % IDENTITY
        U % Coordinate_PS  =  CoordinateUnit
        !-- Local
        U % BaryonMass               =  UNIT % GRAM
        U % NumberDensity            =  UNIT % NUMBER_DENSITY_CGS
        U % MassDensity              =  UNIT % MASS_DENSITY_CGS
        U % EnergyDensity            =  UNIT % ENERGY_DENSITY_CGS
        U % Temperature              =  UNIT % KELVIN
        U % Velocity_U ( 1 )         =  U % Coordinate_PS ( 1 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_CGS
        U % Velocity_U ( 2 )         =  U % Coordinate_PS ( 2 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_CGS
        U % Velocity_U ( 3 )         =  U % Coordinate_PS ( 3 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_CGS
        U % MomentumDensity_D ( 1 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 1 ) ** 2  &
                                        *  U % Velocity_U ( 1 )
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 2 ) ** 2  &
                                        *  U % Velocity_U ( 2 )
        U % MomentumDensity_D ( 3 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 3 ) ** 2  &
                                        *  U % Velocity_U ( 3 )
        !-- Global
        U % Number           =  UNIT % MOLE
        U % Mass             =  UNIT % GRAM
        U % Energy           =  UNIT % ERG
        U % Momentum         =  UNIT % GRAM  *  UNIT % SPEED_CGS
        U % AngularMomentum  =  U % Momentum  *  UNIT % CENTIMETER
      case ( 'ASTROPHYSICS' )
        !-- Spacetime 
        U % Time           =  UNIT % SECOND
        U % Length         =  UNIT % KILOMETER
        U % SqrtDet_M      =  UNIT % IDENTITY
        U % Coordinate_PS  =  CoordinateUnit
        !-- Local
        U % BaryonMass               =  UNIT % MEGA_ELECTRON_VOLT
        U % NumberDensity            =  UNIT % NUMBER_DENSITY_NUCLEAR
        U % MassDensity              =  UNIT % MASS_DENSITY_CGS
        U % EnergyDensity            =  UNIT % ENERGY_DENSITY_NUCLEAR
        U % Temperature              =  UNIT % MEGA_ELECTRON_VOLT
        U % Velocity_U ( 1 )         =  U % Coordinate_PS ( 1 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_OF_LIGHT
        U % Velocity_U ( 2 )         =  U % Coordinate_PS ( 2 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_OF_LIGHT
        U % Velocity_U ( 3 )         =  U % Coordinate_PS ( 3 )  &
                                        /  U % Length  &
                                        *  UNIT % SPEED_OF_LIGHT
        U % MomentumDensity_D ( 1 )  =  UNIT % MASS_DENSITY_NUCLEAR  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 1 ) ** 2  &
                                        *  U % Velocity_U ( 1 )
        U % MomentumDensity_D ( 2 )  =  UNIT % MASS_DENSITY_NUCLEAR  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 2 ) ** 2  &
                                        *  U % Velocity_U ( 2 )
        U % MomentumDensity_D ( 3 )  =  UNIT % MASS_DENSITY_NUCLEAR  &
                                        *  U % Length ** 2  &
                                        /  U % Coordinate_PS ( 3 ) ** 2  &
                                        *  U % Velocity_U ( 3 )
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

  end subroutine Initialize


  subroutine Show_U ( U, IgnorabilityOption )

    class ( Units_F_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call Show ( 'Units -- Spacetime', IgnorabilityOption )
    call Show ( U % Time,          'Time',          IgnorabilityOption )
    call Show ( U % Length,        'Length',        IgnorabilityOption )
    call Show ( U % SqrtDet_M,     'SqrtDet_M',     IgnorabilityOption )
    call Show ( U % Coordinate_PS, 'Coordinate_PS', IgnorabilityOption )

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
    if ( allocated ( U % Coordinate_PS ) ) &
      deallocate ( U % Coordinate_PS )

  end subroutine Finalize


end module Units_F__Form
