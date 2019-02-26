module MarshakWave_Form

  !-- Vaytet et al. 2011

  use GenASiS
  use InteractionsExamples_ASC__Form
  use InteractionsExamples_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( RadiationBoxForm ) :: MarshakWaveForm
  contains
    procedure, private, pass :: &
      Initialize_MW
    generic, public :: &
      Initialize => Initialize_MW
    final :: &
      Finalize
  end type MarshakWaveForm

    private :: &
      InitializeRadiationBox

    ! real ( KDR ), private :: &
    !   AdiabaticIndex, &
    !   SpecificHeatCapacity, &
    !   MeanMolecularWeight, &
    !   MassDensity, &
    !   Temperature, &
    !   TemperatureInner, &
    !   SpecificOpacity, &
    !   SpecificOpacityFloor, &
    !   EnergyScale, &
    !   EnergyMax, &
    !   SoundSpeed

    class ( MarshakWaveForm ), private, pointer :: &
      MarshakWave => null ( )

contains


  subroutine Initialize_MW ( MW, MomentsType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    if ( MW % Type == '' ) &
      MW % Type = 'a MarshakWave'

    MarshakWave => MW

    !-- Parameters

    !-- Initialization

    call InitializeRadiationBox ( MW, MomentsType, Name )

  end subroutine Initialize_MW


  subroutine Finalize ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

  end subroutine Finalize


  subroutine InitializeRadiationBox ( MW, MomentsType, Name )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    integer ( KDI ) :: &
      iD
    real ( KDR ) :: &
      FinishTime
    type ( Character_1D_Form ), dimension ( 3 ) :: &
      BoundaryConditionsFace

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND

    associate ( BCF => BoundaryConditionsFace )
    do iD = 1, 3
      call BCF ( iD ) % Initialize ( [ 'INFLOW', 'INFLOW' ] )     
    end do

    call MW % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'GENERIC' ], &
             MomentsType = MomentsType, &
             Name = Name, &
             BoundaryConditionsFaceOption = BCF, &
             FinishTimeOption = FinishTime )
    end associate !-- BCF

  end subroutine InitializeRadiationBox


end module MarshakWave_Form
