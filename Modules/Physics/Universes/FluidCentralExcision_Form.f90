module FluidCentralExcision_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ApplyGravity_F__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: &
    FluidCentralExcisionForm
      logical ( KDL ) :: &
        Dimensionless = .false.
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FluidCentralExcisionForm

    private :: &
      ApplySources

contains


  subroutine Initialize &
               ( FCE, Name, FluidType, GeometryType, DimensionlessOption, &
                 TimeUnitOption, FinishTimeOption, LimiterParameterOption, &
                 RadiusMaxOption, RadiusMinOption, RadialRatioOption, &
                 CentralMassOption, nCellsPolarOption )

    class ( FluidCentralExcisionForm ), intent ( inout ), target :: &
      FCE
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusMinOption, &
      RadialRatioOption, &
      CentralMassOption, &
      LimiterParameterOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    real ( KDR ) :: &
      RadiusMin, &
      RadiusMax, &
      RadialRatio, &
      CentralMass, &
      FinishTime
    type ( MeasuredValueForm ) :: &
      TimeUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit


    if ( FCE % Type == '' ) &
      FCE % Type = 'a FluidCentralExcision'

    FCE % Dimensionless = .false.
    if ( present ( DimensionlessOption ) ) &
      FCE % Dimensionless = DimensionlessOption


    !-- PositionSpace

    allocate ( Atlas_SC_CE_Form :: FCE % PositionSpace )
    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_CE_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )


    if ( FCE % Dimensionless ) then

      call PS % CreateChart_CE &
             ( RadiusMaxOption = RadiusMaxOption, &
               RadiusMinOption = RadiusMinOption, &
               RadialRatioOption = RadialRatioOption )

    else

      RadiusMax    =  1.0e3_KDR  *  UNIT % KILOMETER
      RadiusMin    =   40.0_KDR  *  UNIT % KILOMETER
      RadialRatio  =  1.0_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax = RadiusMaxOption
      if ( present ( RadiusMinOption ) ) &
        RadiusMin = RadiusMinOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio = RadialRatioOption

      CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

      call PS % CreateChart_CE &
             ( CoordinateUnitOption = CoordinateUnit, &
               RadiusMaxOption = RadiusMax, &
               RadiusMinOption = RadiusMin, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    end if !-- Dimensionless


    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    if ( FCE % Dimensionless ) then
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'CENTRAL_MASS', &
               GravitationalConstantOption = 1.0_KDR, &
               CentralMassOption = 1.0_KDR )
    else
      CentralMass  =  1.4_KDR  *  UNIT % SOLAR_MASS
      if ( present ( CentralMassOption ) ) &
        CentralMass = CentralMassOption
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'CENTRAL_MASS', &
               CentralMassUnitOption = UNIT % SOLAR_MASS, &
               CentralMassOption = CentralMass )
    end if !-- Dimensionless

    call PS % SetGeometry ( GA )


    !-- Fluid

    allocate ( Fluid_ASC_Form :: FCE % Current_ASC )
    select type ( FA => FCE % Current_ASC )
    class is ( Fluid_ASC_Form )

    if ( .not. FCE % Dimensionless ) &
      TimeUnit = UNIT % SECOND
    if ( present ( TimeUnitOption ) ) &
      TimeUnit = TimeUnitOption

    if ( FCE % Dimensionless ) then
      call FA % Initialize ( PS, FluidType )
    else
      call FA % Initialize &
             ( PS, FluidType, &
               Velocity_U_UnitOption &
                 =  CoordinateUnit / TimeUnit, &
               BaryonMassUnitOption &
                 =  UNIT % ATOMIC_MASS_UNIT, &
               NumberDensityUnitOption &
                 =  UNIT % FEMTOMETER ** ( -3 ), &
!                 =  UNIT % MASS_DENSITY_CGS, &
               EnergyDensityUnitOption &
                 =  UNIT % MEGA_ELECTRON_VOLT  &
                    *  UNIT % FEMTOMETER ** ( -3 ), &
!                 =  UNIT % MASS_DENSITY_CGS  *  UNIT % SPEED_OF_LIGHT ** 2, &
               NumberUnitOption &
                 =  UNIT % SOLAR_BARYON_NUMBER, &
!                 =  UNIT % SOLAR_MASS, &
               EnergyUnitOption &
                 =  UNIT % SOLAR_MASS  *  UNIT % SPEED_OF_LIGHT ** 2, &
               MomentumUnitOption &
                 =  UNIT % SOLAR_MASS  *  UNIT % SPEED_OF_LIGHT, &
               AngularMomentumUnitOption &
                 =  UNIT % SOLAR_KERR_PARAMETER, &
!                  =  UNIT % KILOMETER  *  UNIT % SOLAR_MASS  &
!                     *  UNIT % SPEED_OF_LIGHT, &
               TimeUnitOption = TimeUnit, &
               BaryonMassReferenceOption = CONSTANT % ATOMIC_MASS_UNIT, &
               LimiterParameterOption = LimiterParameterOption )
    end if


    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FCE % Step )
    select type ( S => FCE % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources
    end select !-- S


    !-- Template

    FinishTime = 1.0_KDR * TimeUnit
    if ( present ( FinishTimeOption ) ) &
      FinishTime = FinishTimeOption

    call FCE % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, &
             FinishTimeOption = FinishTime )
    call Show ( FCE % Dimensionless, 'Dimensionless' )


    !-- Cleanup

    end select !-- FA
    end select !-- GA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCE )

    type ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FCE

    call FCE % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


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


end module FluidCentralExcision_Form
