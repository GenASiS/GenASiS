module MarshakWave_Form

  use Basics
  use Mathematics
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

    !-- RadiationMoments ( Generic )

    allocate &
      ( RadiationMoments_ASC_Form :: &
          MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    select type ( RMA => MW % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )

    !-- Cleanup

    end select !-- RMA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( MW )
    
    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    call MW % FinalizeTemplate_C_1D ( )

  end subroutine Finalize


end module MarshakWave_Form
