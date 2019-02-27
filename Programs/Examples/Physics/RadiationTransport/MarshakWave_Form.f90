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
      InitializeRadiationBox, &
      InitializeInteractions

    real ( KDR ), private :: &
      BoxLength, &
    !   AdiabaticIndex, &
    !   SpecificHeatCapacity, &
    !   MeanMolecularWeight, &
    !   MassDensity, &
    !   Temperature, &
      TemperatureInner, &
    !   SpecificOpacity, &
    !   SpecificOpacityFloor, &
    !   EnergyScale, &
    !   EnergyMax, &
    !   SoundSpeed
      FinishTime
    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate

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

    associate &
      ( L      => BoxLength, &
        ! Gamma  => AdiabaticIndex, &
        ! C_V    => SpecificHeatCapacity, &
        ! N_0    => MassDensity, &
        ! T_0    => Temperature, &
        T_I    => TemperatureInner )!, &
        ! Kappa  => SpecificOpacity, &
        ! Kappa_Min => SpecificOpacityFloor )

    L      =  20.0_KDR    *  UNIT % CENTIMETER
    ! Gamma  =  1.4_KDR
    ! C_V    =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    ! N_0    =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    ! T_0    =  3.0e2_KDR   *  UNIT % KELVIN
    T_I    =  1.0e3_KDR   *  UNIT % KELVIN
    ! Kappa  =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    ! Kappa_Min    =  10.0_KDR  *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    ! EnergyScale  =  T_I

    call PROGRAM_HEADER % GetParameter ( L,     'BoxLength' )
    ! call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    ! call PROGRAM_HEADER % GetParameter ( C_V,   'SpecificHeatCapacity' )
    ! call PROGRAM_HEADER % GetParameter ( N_0,   'MassDensity' )
    ! call PROGRAM_HEADER % GetParameter ( T_0,   'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,   'TemperatureInner' )
    ! call PROGRAM_HEADER % GetParameter ( Kappa, 'SpecificOpacity' )
    ! call PROGRAM_HEADER % GetParameter ( Kappa, 'SpecificOpacity' )
    ! call PROGRAM_HEADER % GetParameter ( Kappa_Min,   'SpecificOpacityFloor' )
    ! call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND


    !-- Initialization

    call InitializeRadiationBox ( MW, MomentsType, Name )


    !-- Cleanup

    end associate !-- L, etc.

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
      EnergyScale
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit_PS, &
      CoordinateUnit_MS, &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( Character_1D_Form ), dimension ( 3 ) :: &
      BoundaryConditionsFace


    !-- Position space parameters

    associate ( BCF => BoundaryConditionsFace )
    do iD = 1, 3
      call BCF ( iD ) % Initialize ( [ 'INFLOW', 'INFLOW' ] )     
    end do

    MinCoordinate  =  0.0_KDR
    MaxCoordinate  =  BoxLength


    !-- Momentum space parameters

    EnergyScale = TemperatureInner


    !-- Units

    CoordinateUnit_PS  =  UNIT % CENTIMETER

    CoordinateUnit_MS = UNIT % IDENTITY
    CoordinateUnit_MS ( 1 ) = UNIT % ELECTRON_VOLT

    TimeUnit = UNIT % SECOND

    BaryonMassUnit     =  UNIT % GRAM
    NumberDensityUnit  =  UNIT % MOLE  *  UNIT % CENTIMETER ** (-3)
    EnergyDensityUnit  =  UNIT % ERG   *  UNIT % CENTIMETER ** (-3)
    TemperatureUnit    =  UNIT % KELVIN

    Velocity_U_Unit  =  UNIT % CENTIMETER  /  UNIT % SECOND 

    MomentumDensity_U_Unit &
      =  UNIT % GRAM  *  UNIT % CENTIMETER ** (-2)  /  UNIT % SECOND
    MomentumDensity_D_Unit &
      =  UNIT % GRAM  *  UNIT % CENTIMETER ** (-2)  /  UNIT % SECOND


    !-- Initialization

    call MW % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'GENERIC' ], &
             MomentsType = MomentsType, &
             Name = Name, &
             BoundaryConditionsFaceOption = BCF, &
             CoordinateUnit_PS_Option = CoordinateUnit_PS, &
             CoordinateUnit_MS_Option = CoordinateUnit_MS, &
             Velocity_U_UnitOption = Velocity_U_Unit, &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             TimeUnitOption = TimeUnit, &
             BaryonMassUnitOption = BaryonMassUnit, &
             NumberDensityUnitOption = NumberDensityUnit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             TemperatureUnitOption = TemperatureUnit, &
             FinishTimeOption = FinishTime, &
             EnergyScaleOption = EnergyScale, &
             BaryonMassReferenceOption = CONSTANT % ATOMIC_MASS_UNIT )
    end associate !-- BCF

  end subroutine InitializeRadiationBox


  subroutine InitializeInteractions ( MW )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    ! select type ( I => T % Integrator )
    ! class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

    !   select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
    !   class is ( RadiationMoments_ASC_Form )

    !   select type ( FA => I % Current_ASC )
    !   class is ( Fluid_ASC_Form )

    !   select type ( PS => I % PositionSpace )
    !   class is ( Atlas_SC_Form )

    !   allocate ( InteractionsExamples_ASC_Form :: T % Interactions_ASC )
    !   select type ( IA => T % Interactions_ASC )
    !   class is ( InteractionsExamples_ASC_Form )
    !   call IA % Initialize &
    !          ( PS, InteractionsType = 'CONSTANT', MomentsType = 'GREY' )
    !   call IA % Set ( FA, OpacityAbsorption = OpacityAbsorption )
    !   call RMA % SetInteractions ( IA )
    !   end select !-- IA

    !   end select !-- PS
    !   end select !-- FA
    !   end select !-- RMA

    ! class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

    !   select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    !   class is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    !   select type ( FA => I % Current_ASC )
    !   class is ( Fluid_ASC_Form )

    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )

    !   allocate ( InteractionsExamples_BSLL_ASC_CSLD_Form :: &
    !                T % Interactions_BSLL_ASC_CSLD )
    !   select type ( IB => T % Interactions_BSLL_ASC_CSLD )
    !   class is ( InteractionsExamples_BSLL_ASC_CSLD_Form )
    !   call IB % Initialize ( MS, InteractionsType = 'CONSTANT' )
    !   call IB % Set ( FA, OpacityAbsorption = OpacityAbsorption )
    !   call RMB % SetInteractions ( IB )
    !   end select !-- IB

    !   end select !-- PS
    !   end select !-- FA
    !   end select !-- RMA

    ! end select !-- Integrator

  end subroutine InitializeInteractions


end module MarshakWave_Form
