module MarshakWave_Form

  use Basics
  use Mathematics
  use Fluid_P_NR__Form
  use Fluid_ASC__Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use Interactions_BE_G__Form
  use Interactions_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_1D_Template ) :: MarshakWaveForm
    integer ( KDI ) :: &
      RADIATION = 1, &
      FLUID     = 2         
    type ( Interactions_ASC_Form ), allocatable :: &
      Interactions_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type MarshakWaveForm

    private :: &
      SetRadiation, &
      SetFluid, &
      SetInteractions

    real ( KDR ), private :: &
      AdiabaticIndex, &
      MassDensity, &
      Temperature, &
      TemperatureInner, &
      SpecificOpacity
    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate

contains


  subroutine Initialize ( MW, Name )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    real ( KDR ) :: &
      FinishTime, &
TimeStep
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      MassDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit, &
      MassUnit, &
      EnergyUnit, &
      MomentumUnit, &
      AngularMomentumUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit, &
      VelocityUnit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit

    if ( MW % Type == '' ) &
      MW % Type = 'a MarshakWave' 

    associate &
      ( Gamma => AdiabaticIndex, &
        N_0   => MassDensity, &
        T_0   => Temperature, &
        T_I   => TemperatureInner, &
        Kappa => SpecificOpacity )

    Gamma  =  1.4_KDR
    N_0    =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_0    =  3.0e2_KDR   *  UNIT % KELVIN
    T_I    =  1.0e3_KDR   *  UNIT % KELVIN
    Kappa  =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( N_0,   'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,   'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,   'TemperatureInner' )
    call PROGRAM_HEADER % GetParameter ( Kappa, 'SpecificOpacity' )

    call Show ( 'MarshakWave parameters' )
    call Show ( Gamma, 'Gamma' )
    call Show ( N_0, UNIT % MASS_DENSITY_CGS, 'N_0' )
    call Show ( T_0, UNIT % KELVIN, 'T_0' )
    call Show ( T_I, UNIT % KELVIN, 'T_I' )
    call Show ( Kappa, UNIT % CENTIMETER ** 2 / UNIT % GRAM, 'Kappa' )

    end associate !-- Gamma, etc.

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: MW % PositionSpace )
    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    do iD = 1, PS % nDimensions
      call PS % SetBoundaryConditionsFace &
             ( [ 'INFLOW', 'INFLOW' ], iDimension = iD )
    end do !-- iD

    CoordinateUnit  =  UNIT % CENTIMETER

    MinCoordinate  =   0.0_KDR  *  UNIT % CENTIMETER
    MaxCoordinate  =  20.0_KDR  *  UNIT % CENTIMETER

    call PS % CreateChart &
           ( CoordinateUnitOption = CoordinateUnit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate )

    !-- Geometry of PositionSpace

    allocate ( MW % Geometry_ASC )
    associate ( GA => MW % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Prepare for Currents

    MW % N_CURRENTS = 2
    allocate ( MW % Current_ASC_1D ( MW % N_CURRENTS ) )
    allocate ( MW % TimeStepLabel ( MW % N_CURRENTS ) )
    MW % TimeStepLabel ( MW % RADIATION ) = 'Radiation'
    MW % TimeStepLabel ( MW % FLUID )     = 'Fluid'

    TimeUnit = UNIT % SECOND

    VelocityUnit ( 1 )     =  CoordinateUnit ( 1 ) / TimeUnit 
    VelocityUnit ( 2 )     =  CoordinateUnit ( 2 ) / TimeUnit
    VelocityUnit ( 3 )     =  CoordinateUnit ( 3 ) / TimeUnit
    MassDensityUnit        =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit      =  UNIT % MASS_DENSITY_CGS  &
                              *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumDensity_U_Unit = EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
    MomentumDensity_D_Unit = EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
    TemperatureUnit        =  UNIT % KELVIN

    MassUnit  =  MassDensityUnit  *  CoordinateUnit ( 1 )
    if ( PS % nDimensions > 1 ) &
      MassUnit  =  MassUnit  *  CoordinateUnit ( 2 )
    if ( PS % nDimensions > 2 ) &
      MassUnit  =  MassUnit  *  CoordinateUnit ( 3 )

    EnergyUnit           =  MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumUnit         =  MassUnit  *  UNIT % SPEED_OF_LIGHT
    AngularMomentumUnit  =  CoordinateUnit ( 1 ) &
                            *  MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    !-- RadiationMoments

    allocate &
      ( RadiationMoments_ASC_Form :: &
          MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    select type ( RA => MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA % Initialize &
           ( PS, 'GENERIC', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = EnergyDensityUnit )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: MW % Current_ASC_1D ( MW % FLUID ) % Element )
    select type ( FA => MW % Current_ASC_1D ( MW % FLUID ) % Element )
    class is ( Fluid_ASC_Form )

    call FA % Initialize &
           ( PS, 'NON_RELATIVISTIC', VelocityUnitOption = VelocityUnit, &
             MassDensityUnitOption = MassDensityUnit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             TemperatureUnitOption = TemperatureUnit, &
             MassUnitOption = MassUnit, EnergyUnitOption = EnergyUnit, &
             MomentumUnitOption = MomentumUnit, &
             AngularMomentumUnitOption = AngularMomentumUnit )

    !-- Interactions

    allocate ( MW % Interactions_ASC )
    associate ( IA => MW % Interactions_ASC )
    call IA % Initialize &
           ( PS, 'BOSE_EINSTEIN_GREY', &
             LengthUnitOption = CoordinateUnit ( 1 ), &
             EnergyDensityUnitOption = EnergyDensityUnit )
    call RA % SetInteractions ( IA )

    !-- Initial conditions

    call SetRadiation ( MW )
    call SetFluid ( MW )
    call SetInteractions ( MW )

    !-- Initialize template

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND
 
    call MW % InitializeTemplate_C &
           ( Name, TimeUnitOption = TimeUnit, FinishTimeOption = 0.0_KDR )

    call MW % ComputeTimeStep ( TimeStep )

    !-- Cleanup

    end associate !-- IA
    end select !-- FA
    end select !-- RA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( MW )
    
    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    if ( allocated ( MW % Interactions_ASC) ) &
      deallocate ( MW % Interactions_ASC )

    call MW % FinalizeTemplate_C_1D ( )

  end subroutine Finalize


  subroutine SetRadiation ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( RadiationMomentsForm ), pointer :: &
      R

    associate &
      ( T_0 => Temperature, &
        T_I => TemperatureInner )

    select type ( RA => MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    class is ( RadiationMoments_ASC_Form )
    R => RA % RadiationMoments ( )

    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (   J => R % Value ( :, R % COMOVING_ENERGY_DENSITY ), &
        H_1 => R % Value ( :, R % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
        H_2 => R % Value ( :, R % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
        H_3 => R % Value ( :, R % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
          X => G % Value ( :, G % CENTER ( 1 ) ), &
          Y => G % Value ( :, G % CENTER ( 2 ) ), &
          Z => G % Value ( :, G % CENTER ( 3 ) ), &
          a => CONSTANT % RADIATION )

    J  =  a  *  T_0 ** 4

    H_1  =  0.0_KDR
    H_2  =  0.0_KDR
    H_3  =  0.0_KDR

    where ( X < MinCoordinate ( 1 ) .or. Y < MinCoordinate ( 2 ) &
            .or. Z < MinCoordinate ( 3 ) )
      J  =  a  *  T_I ** 4
    end where

    call R % ComputeFromPrimitive ( G )

    end associate !-- J, etc.
    end select !-- PS
    end select !-- RA
    end associate !-- T_0
    nullify ( R, G )

  end subroutine SetRadiation


  subroutine SetFluid ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_NR_Form ), pointer :: &
      F

    associate &
      ( Gamma => AdiabaticIndex, &
        N_0   => MassDensity, &
        T_0   => Temperature )

    select type ( FA => MW % Current_ASC_1D ( MW % FLUID ) % Element )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_NR ( )

    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call F % SetAdiabaticIndex ( Gamma )
    call F % SetFiducialBaryonDensity ( N_0 )
    call F % SetFiducialTemperature ( T_0 )

    associate &
      (   N => F % Value ( :, F % COMOVING_DENSITY ), &
        V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
        V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
          T => F % Value ( :, F % TEMPERATURE ) )

    N  =  N_0
    T  =  T_0

    V_1  =  0.0_KDR
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR

    call F % ComputeFromPrimitiveTemperature ( G )

    end associate !-- N, etc.
    end select !-- PS
    end select !-- FA
    end associate !-- Gamma, etc.
    nullify ( F, G )

  end subroutine SetFluid


  subroutine SetInteractions ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    class ( Fluid_P_NR_Form ), pointer :: &
      F
    class ( Interactions_BE_G_Form ), pointer :: &
      I

    associate ( IA => MW % Interactions_ASC )

    select type ( FA => MW % Current_ASC_1D ( MW % FLUID ) % Element )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_NR ( )

    I => IA % Interactions_BE_G ( )
    call I % Set ( F, SpecificOpacity )
    call I % Compute ( )

    end select !-- FA
    end associate !-- IA
    nullify ( F, I )

  end subroutine SetInteractions


end module MarshakWave_Form
