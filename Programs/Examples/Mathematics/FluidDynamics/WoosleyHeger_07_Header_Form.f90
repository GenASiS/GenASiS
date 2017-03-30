module WoosleyHeger_07_Header_Form

  use Basics
  use Mathematics
  use Fluid_ASC__Form

  implicit none
  private

  type, public :: WoosleyHeger_07_HeaderForm
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
      VelocityUnit
    class ( IntegratorTemplate ), pointer :: &
      Integrator => null ( )
    class ( Fluid_ASC_Form ), pointer :: &
      Fluid_ASC => null ( )
  contains
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeFluid
  end type WoosleyHeger_07_HeaderForm

contains


  subroutine InitializePositionSpace ( WHH, I )

    class ( WoosleyHeger_07_HeaderForm ), intent ( inout ) :: &
      WHH
    class ( IntegratorTemplate ), intent ( in ), target :: &
      I

    integer ( KDI ) :: &
      nCellsCore, &
      nCellsRadius, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    WHH % Integrator => I

    allocate ( Atlas_SC_Form :: WHH % Integrator % PositionSpace )
    select type ( PS => WHH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    CoordinateSystem = 'SPHERICAL'
    WHH % CoordinateUnit   = [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

    call PS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
    if ( PS % nDimensions > 1 ) &
      call PS % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    if ( PS % nDimensions > 2 ) &
      call PS % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 3 )

    RadiusMax = 1.0e4_KDR * UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'RadiusMax' )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [   0.0_KDR, 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    RadiusCore = 40.0_KDR  *  UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( RadiusCore, 'RadiusCore' )

    nCellsCore = 32  !-- Number of central cells with equal spacing
    call PROGRAM_HEADER % GetParameter ( nCellsCore, 'nCellsCore' )

    call Show ( 'Mesh core parameters' )
    call Show ( RadiusCore, UNIT % KILOMETER, 'RadiusCore' )
    call Show ( nCellsCore, 'nCellsCore' )
    call Show ( RadiusCore / nCellsCore, UNIT % KILOMETER, 'CellWidthCore' )

    nCellsRadius = 6.5 * nCellsCore  !-- Aiming for roughly 10,000 km
    call PROGRAM_HEADER % GetParameter ( nCellsRadius, 'nCellsRadius' )

    nCellsPolar     = 3 * nCellsCore
    nCellsAzimuthal = 2 * nCellsPolar

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  RadiusCore

    nCells = [ nCellsRadius, 1, 1 ]
    if ( PS % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( PS % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    call PS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = WHH % CoordinateUnit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nEqualOption = nCellsCore )

    !-- Geometry of PositionSpace

    allocate ( WHH % Integrator % Geometry_ASC )
    associate ( GA => WHH % Integrator % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    WHH % TimeUnit = UNIT % SECOND

    end select !-- PS

  end subroutine InitializePositionSpace


  subroutine InitializeFluid ( WHH, FA )

    class ( WoosleyHeger_07_HeaderForm ), intent ( inout ) :: &
      WHH
    type ( Fluid_ASC_Form ), intent ( inout ), target :: &
      FA

    WHH % Fluid_ASC => FA

    WHH % VelocityUnit       =  WHH % CoordinateUnit  /  WHH % TimeUnit 
    WHH % MassDensityUnit    =  UNIT % MASS_DENSITY_CGS
    WHH % EnergyDensityUnit  =  UNIT % MASS_DENSITY_CGS  &
                                *  UNIT % SPEED_OF_LIGHT ** 2
    WHH % TemperatureUnit    =  UNIT % MEV

    WHH % MassUnit             =  UNIT % SOLAR_MASS
    WHH % EnergyUnit           =  WHH % MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    WHH % MomentumUnit         =  WHH % MassUnit  *  UNIT % SPEED_OF_LIGHT
    WHH % AngularMomentumUnit  =  WHH % CoordinateUnit ( 1 ) &
                                  *  WHH % MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    select type ( PS => WHH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    call FA % Initialize &
           ( PS, 'MEAN_HEAVY_NUCLEUS', &
             VelocityUnitOption        = WHH % VelocityUnit, &
             MassDensityUnitOption     = WHH % MassDensityUnit, &
             EnergyDensityUnitOption   = WHH % EnergyDensityUnit, &
             TemperatureUnitOption     = WHH % TemperatureUnit, &
             MassUnitOption            = WHH % MassUnit, &
             EnergyUnitOption          = WHH % EnergyUnit, &
             MomentumUnitOption        = WHH % MomentumUnit, &
             AngularMomentumUnitOption = WHH % AngularMomentumUnit, &
             ShockThresholdOption = 1.0_KDR )

    end select !-- PS

  end subroutine InitializeFluid


end module WoosleyHeger_07_Header_Form
