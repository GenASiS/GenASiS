module MarshakWave_Form

  use Basics
  use Mathematics
  use Fluid_ASC__Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_1D_Template ) :: MarshakWaveForm
    integer ( KDI ) :: &
      RADIATION = 1, &
      FLUID     = 2         
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type MarshakWaveForm

contains


  subroutine Initialize ( MW, Name )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      MassDensityUnit, &
      EnergyDensityUnit, &
      MassUnit, &
      EnergyUnit, &
      MomentumUnit, &
      AngularMomentumUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit, &
      VelocityUnit

    if ( MW % Type == '' ) &
      MW % Type = 'a MarshakWave' 

    !-- PositionSpace

    CoordinateUnit  =  UNIT % CENTIMETER

    MinCoordinate  =   0.0_KDR  *  UNIT % CENTIMETER
    MaxCoordinate  =  20.0_KDR  *  UNIT % CENTIMETER

    allocate ( Atlas_SC_Form :: MW % PositionSpace )
    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
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

    !-- RadiationMoments

    allocate &
      ( RadiationMoments_ASC_Form :: &
          MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    select type ( RMA => MW % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )

    !-- Fluid

    TimeUnit = UNIT % SECOND

    VelocityUnit ( 1 ) =  CoordinateUnit ( 1 ) / TimeUnit 
    VelocityUnit ( 2 ) =  CoordinateUnit ( 2 ) / TimeUnit
    VelocityUnit ( 3 ) =  CoordinateUnit ( 3 ) / TimeUnit
    MassDensityUnit    =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit  =  UNIT % MASS_DENSITY_CGS  &
                          *  UNIT % SPEED_OF_LIGHT ** 2

    MassUnit             =  UNIT % GRAM
    EnergyUnit           =  MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumUnit         =  MassUnit  *  UNIT % SPEED_OF_LIGHT
    AngularMomentumUnit  =  CoordinateUnit ( 1 ) &
                            *  MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    allocate ( Fluid_ASC_Form :: MW % Current_ASC_1D ( MW % FLUID ) % Element )
    select type ( FA => MW % Current_ASC )
    class is ( Fluid_ASC_Form )

    call FA % Initialize &
           ( PS, 'NON_RELATIVISTIC', VelocityUnitOption = VelocityUnit, &
             MassDensityUnitOption = MassDensityUnit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             MassUnitOption = MassUnit, EnergyUnitOption = EnergyUnit, &
             MomentumUnitOption = MomentumUnit, &
             AngularMomentumUnitOption = AngularMomentumUnit )

    !-- Cleanup

    end select !-- FA
    end select !-- RMA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( MW )
    
    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    call MW % FinalizeTemplate_C_1D ( )

  end subroutine Finalize


end module MarshakWave_Form
