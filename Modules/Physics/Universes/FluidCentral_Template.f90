module FluidCentral_Template

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ApplyGravity_F__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    FluidCentralTemplate
      type ( MeasuredValueForm ), dimension ( 3 ) :: &
        CoordinateUnit
      logical ( KDL ) :: &
        Dimensionless
  contains
    procedure, public, pass :: &
      InitializeTemplate_FC
    procedure ( IPS ), private, pass, deferred :: &
      InitializePositionSpace
    procedure ( IG ), private, pass, deferred :: &
      InitializeGeometry
  end type FluidCentralTemplate

  abstract interface

    subroutine IPS ( FC, RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
                     RadialRatioOption, nCellsCoreOption, nCellsPolarOption )
      use Basics
      import FluidCentralTemplate
      class ( FluidCentralTemplate ), intent ( inout ) :: &
        FC
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusMinOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsCoreOption, &
        nCellsPolarOption
    end subroutine IPS

    subroutine IG ( FC, GA, PS, GeometryType, CentralMassOption )
      use Basics
      use Mathematics
      use Spaces
      import FluidCentralTemplate
      class ( FluidCentralTemplate ), intent ( inout ) :: &
        FC
      type ( Geometry_ASC_Form ), intent ( inout ) :: &
        GA
      class ( Atlas_SC_Form ), intent ( in ) :: &
        PS
      character ( * ), intent ( in )  :: &
        GeometryType
      real ( KDR ), intent ( in ), optional :: &
        CentralMassOption
    end subroutine IG

  end interface

    private :: &
      ApplySources

contains


  subroutine InitializeTemplate_FC &
               ( FC, Name, FluidType, GeometryType, DimensionlessOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, &
                 RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
                 RadialRatioOption, CentralMassOption, nWriteOption, &
                 nCellsCoreOption, nCellsPolarOption )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      LimiterParameterOption, &
      ShockThresholdOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusMinOption, &
      RadialRatioOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption, &
      nCellsCoreOption, &
      nCellsPolarOption

    real ( KDR ) :: &
      FinishTime
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_D_Unit


    if ( FC % Type == '' ) &
      FC % Type = 'a FluidCentral'

    FC % Dimensionless = .false.
    if ( present ( DimensionlessOption ) ) &
      FC % Dimensionless = DimensionlessOption

    if ( .not. FC % Dimensionless ) then
      FC % CoordinateUnit  &
        =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]
      TimeUnit &
        = UNIT % SECOND
    end if

    if ( present ( TimeUnitOption ) ) &
      TimeUnit = TimeUnitOption


    !-- PositionSpace

    call FC % InitializePositionSpace &
           ( RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
             RadialRatioOption, nCellsCoreOption, nCellsPolarOption )

    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    call FC % InitializeGeometry ( GA, PS, GeometryType, CentralMassOption )

    call PS % SetGeometry ( GA )


    !-- Fluid

    allocate ( Fluid_ASC_Form :: FC % Current_ASC )
    select type ( FA => FC % Current_ASC )
    class is ( Fluid_ASC_Form )

    if ( FC % Dimensionless ) then
      call FA % Initialize ( PS, FluidType )
    else
      BaryonMassUnit                =  UNIT % ATOMIC_MASS_UNIT
      NumberDensityUnit             =  UNIT % NUMBER_DENSITY_NUCLEAR
      EnergyDensityUnit             =  UNIT % ENERGY_DENSITY_NUCLEAR
      Velocity_U_Unit               =  FC % CoordinateUnit  /  TimeUnit
      MomentumDensity_D_Unit        =  BaryonMassUnit * NumberDensityUnit &
                                       * Velocity_U_Unit
      MomentumDensity_D_Unit ( 2 )  =  MomentumDensity_D_Unit ( 2 ) &
                                       *  FC % CoordinateUnit ( 1 ) ** 2
      MomentumDensity_D_Unit ( 3 )  =  MomentumDensity_D_Unit ( 3 ) &
                                       *  FC % CoordinateUnit ( 1 ) ** 2
      call FA % Initialize &
             ( PS, FluidType, &
               Velocity_U_UnitOption         =  Velocity_U_Unit, &
               MomentumDensity_D_UnitOption  =  MomentumDensity_D_Unit, &
               BaryonMassUnitOption          =  BaryonMassUnit, &
               NumberDensityUnitOption       =  NumberDensityUnit , &
               EnergyDensityUnitOption       =  EnergyDensityUnit, &
               NumberUnitOption              =  UNIT % SOLAR_BARYON_NUMBER, &
               EnergyUnitOption              =  UNIT % ENERGY_SOLAR_MASS, &
               MomentumUnitOption            =  UNIT % MOMENTUM_SOLAR_MASS, &
               AngularMomentumUnitOption     =  UNIT % SOLAR_KERR_PARAMETER, &
               TimeUnitOption                =  TimeUnit, &
               BaryonMassReferenceOption     =  CONSTANT % ATOMIC_MASS_UNIT, &
               LimiterParameterOption        =  LimiterParameterOption, &
               ShockThresholdOption          =  ShockThresholdOption )
    end if


    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FC % Step )
    select type ( S => FC % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources
    end select !-- S


    !-- Template

    FinishTime = 1.0_KDR * TimeUnit
    if ( present ( FinishTimeOption ) ) &
      FinishTime = FinishTimeOption

    call FC % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, &
             FinishTimeOption = FinishTime, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

    call Show ( FC % Dimensionless, 'Dimensionless' )


    !-- Cleanup

    end select !-- FA
    end select !-- GA
    end select !-- PS

  end subroutine InitializeTemplate_FC


  subroutine ApplySources &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call ApplyCurvilinear_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    call ApplyGravity_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

  end subroutine ApplySources


end module FluidCentral_Template
