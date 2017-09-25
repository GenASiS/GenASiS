module MarshakWave_Form

  !-- Vaytet et al. 2011

  use Basics
  use Mathematics
  use Fluid_P_NR__Form
  use Sources_F__Form
  use Fluid_ASC__Form
  use SetPlanckSpectrum_Command
  use PhotonMoments_S__Form
  use ApplyRelaxation_RM__Command
  use RadiationMoments_BSLL_ASC_CSLD__Form
  use Interactions_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( Integrator_C_MS_C_PS_Template ) :: MarshakWaveForm
    type ( Interactions_BSLL_ASC_CSLD_Form ), allocatable :: &
      Interactions_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
  end type MarshakWaveForm

    private :: &
      SetRadiation, &
      SetFluid!, &
      ! SetInteractions, &
      ! ApplySources_Fluid

      ! private :: &
      !   ComputeFluidSource_G_S_Radiation_Kernel, &
      !   ApplySources_Fluid_Kernel

    real ( KDR ), private :: &
      AdiabaticIndex, &
      SpecificHeatCapacity, &
      MeanMolecularWeight, &
      MassDensity, &
       Temperature, &
       TemperatureInner, &
    !   SpecificOpacity, &
      EnergyScale, &
      EnergyMax, &
      SoundSpeed
    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDF ) :: &
      InteractionsType
    class ( Fluid_P_NR_Form ), pointer :: &
      Fluid => null ( )
    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), pointer :: &
      RadiationBundle => null ( )
    ! class ( InteractionsTemplate ), pointer :: &
    !   Interactions => null ( )

contains


  subroutine Initialize ( MW, Name )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      nEnergyCells
    real ( KDR ) :: &
      BoxLength, &
    !   MeanFreePath, &
    !   OpticalDepth, &
      DynamicalTime, &
    !   DiffusionTime, &
      FinishTime
    real ( KDR ), dimension ( 3 ) :: &
      Scale
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
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    associate &
      ( L      => BoxLength, &
        Gamma  => AdiabaticIndex, &
        C_V    => SpecificHeatCapacity, &
        N_0    => MassDensity, &
        T_0    => Temperature, &
        T_I    => TemperatureInner &!, &
    !     Kappa  => SpecificOpacity )
    )

    L      =  20.0_KDR    *  UNIT % CENTIMETER
    Gamma  =  1.4_KDR
    C_V    =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    N_0    =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_0    =  3.0e2_KDR   *  UNIT % KELVIN
    T_I    =  1.0e3_KDR   *  UNIT % KELVIN
    ! Kappa  =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    EnergyScale  =  sqrt ( T_0 * T_I )

    call PROGRAM_HEADER % GetParameter ( L,     'BoxLength' )
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,   'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( N_0,   'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,   'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,   'TemperatureInner' )
    ! call PROGRAM_HEADER % GetParameter ( Kappa, 'SpecificOpacity' )
    call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

    InteractionsType = 'MARSHAK_WAVE_VAYTET_1_SPECTRAL'
    call PROGRAM_HEADER % GetParameter ( InteractionsType, 'InteractionsType' )

    select case ( trim ( InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_2_SPECTRAL', 'MARSHAK_WAVE_VAYTET_3_SPECTRAL' )
      EnergyMax  =  0.620_KDR * UNIT % ELECTRON_VOLT
      call PROGRAM_HEADER % GetParameter ( EnergyMax, 'EnergyMax' )
    end select !-- InteractionsType

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: MW % PositionSpace )
    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    do iD = 1, PS % nDimensions
      call PS % SetBoundaryConditionsFace &
             ( [ 'INFLOW', 'INFLOW' ], iDimension = iD )
    end do !-- iD

    CoordinateUnit  =  UNIT % CENTIMETER

    MinCoordinate  =  0.0_KDR
    MaxCoordinate  =  BoxLength

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

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: MW % MomentumSpace )
    select type ( MS => MW % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, Name )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = EnergyScale

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % ELECTRON_VOLT

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    nEnergyCells = 20
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Prepare for Currents

    allocate ( MW % TimeStepLabel ( 2 ) )
    MW % TimeStepLabel ( 1 )  =  'Radiation'
    MW % TimeStepLabel ( 2 )  =  'Fluid'

    TimeUnit = UNIT % SECOND

    VelocityUnit ( 1 )     =  CoordinateUnit ( 1 ) / TimeUnit 
    VelocityUnit ( 2 )     =  CoordinateUnit ( 2 ) / TimeUnit
    VelocityUnit ( 3 )     =  CoordinateUnit ( 3 ) / TimeUnit
    MassDensityUnit        =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit      =  UNIT % MASS_DENSITY_CGS  &
                              *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumDensity_U_Unit =  EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
    MomentumDensity_D_Unit =  EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
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
    
    !-- Radiation

    allocate ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
               MW % Current_BSLL_ASC_CSLD )
    select type ( RMB => MW % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )
    call RMB % Initialize &
           ( MS, 'PHOTONS_SPECTRAL', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = EnergyDensityUnit )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: MW % Current_ASC )
    select type ( FA => MW % Current_ASC )
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

    allocate ( MW % Interactions_BSLL_ASC_CSLD )
    associate ( IB => MW % Interactions_BSLL_ASC_CSLD )
    call IB % Initialize &
           ( MS, InteractionsType, &
             LengthUnitOption = CoordinateUnit ( 1 ), &
             EnergyDensityUnitOption = EnergyDensityUnit )
    call RMB % SetInteractions ( IB )
    end associate !-- IB

    !-- Step

    allocate ( Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form :: MW % Step )
    select type ( S => MW % Step )
    class is ( Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form )
    call S % Initialize ( RMB, FA, Name )
    S % ApplyRelaxation_F % Pointer => ApplyRelaxation_RM
    end select !-- S

    !-- Initial conditions

    call SetRadiation ( MW )
    call SetFluid ( MW )
    ! call SetInteractions ( MW )

    associate &
      ( Mu     => MeanMolecularWeight, & 
        c_s    => SoundSpeed, &
        t_Dyn  => DynamicalTime, &
    !     Lambda => MeanFreePath, &
    !     Tau    => OpticalDepth, &
    !     t_Diff => DiffusionTime, &
        c      => CONSTANT % SPEED_OF_LIGHT &
      )

    t_Dyn   =  L / c_s

    ! Lambda  =  1.0_KDR / ( N_0 * Kappa )
    ! Tau     =  L / Lambda
    ! t_Diff  =  L * Tau / c

    call Show ( 'MarshakWave parameters' )
    call Show ( Gamma, 'Gamma' )
    call Show ( C_V, UNIT % ERG / UNIT % KELVIN / UNIT % GRAM, 'C_V' )
    call Show ( Mu, 'Mu' )
    call Show ( N_0, UNIT % MASS_DENSITY_CGS, 'N_0' )
    call Show ( T_0, UNIT % KELVIN, 'T_0' )
    call Show ( c_s, UNIT % SPEED_OF_LIGHT, 'c_s' )
    call Show ( t_Dyn, UNIT % SECOND, 't_Dyn' )
    call Show ( T_I, UNIT % KELVIN, 'T_I' )
    ! call Show ( Kappa, UNIT % CENTIMETER ** 2 / UNIT % GRAM, 'Kappa' )
    ! call Show ( Lambda, UNIT % CENTIMETER, 'Lambda' )
    ! call Show ( Tau, UNIT % IDENTITY, 'Tau' )
    ! call Show ( t_Diff, UNIT % SECOND, 't_Diff' )
    ! call Show ( InteractionsType, 'InteractionsType' )

    ! select case ( trim ( InteractionsType ) )
    ! case ( 'MARSHAK_WAVE_VAYTET_2_GREY', 'MARSHAK_WAVE_VAYTET_3_GREY' )
    !   call Show ( EnergyMax, UNIT % ELECTRON_VOLT, 'EnergyMax' )
    ! end select !-- InteractionsType

    end associate !-- Mu, etc.

    !-- Initialize template

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND

    call MW % InitializeTemplate_C_MS_C_PS &
           ( Name, TimeUnitOption = TimeUnit, FinishTimeOption = FinishTime )

    !-- Cleanup

    end select !-- FA
    end select !-- RMB
    end select !-- MS
    end select !-- PS
    end associate !-- L, etc.

  end subroutine Initialize


  impure elemental subroutine Finalize ( MW )
    
    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    if ( allocated ( MW % Interactions_BSLL_ASC_CSLD ) ) &
      deallocate ( MW % Interactions_BSLL_ASC_CSLD )

    call MW % FinalizeTemplate_C_MS_C_PS ( )

  end subroutine Finalize


  subroutine SetRadiation ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( PhotonMoments_S_Form ), pointer :: &
      RM

    select type ( MS => MW % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    G => MS % Base_CSLD % Geometry ( )

    select type ( RMB => MW % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    do iF = 1, MS % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )
      RM => RMB % PhotonMoments_S ( iF )
      associate &
        ( J    =>  RM % Value ( :, RM % COMOVING_ENERGY ), &
          H_1  =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
          H_2  =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
          H_3  =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
          T_0  =>  Temperature, &
          T_I  =>  TemperatureInner, &
          E    =>  RMB % Energy, &
          X    =>  G % Value ( iBC, G % CENTER ( 1 ) ), &
          Y    =>  G % Value ( iBC, G % CENTER ( 2 ) ), &
          Z    =>  G % Value ( iBC, G % CENTER ( 3 ) ) )

      if ( X < MinCoordinate ( 1 ) .or. Y < MinCoordinate ( 2 ) &
            .or. Z < MinCoordinate ( 3 ) ) &
      then
        call SetPlanckSpectrum ( E, T_I, J )
      else
        call SetPlanckSpectrum ( E, T_0, J )
      end if

      H_1  =  0.0_KDR
      H_2  =  0.0_KDR
      H_3  =  0.0_KDR

      call RM % ComputeFromPrimitive ( iBC, G )

      end associate !-- J, etc.
      end associate !-- iBC
    end do !-- iF

    call RMB % LoadSections ( )

    !-- Module variable for accessibility
    RadiationBundle => RMB

    end select !-- RMB
    end select !-- MS

    nullify ( G, RM )

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
        C_V   => SpecificHeatCapacity, &
        N_0   => MassDensity, &
        T_0   => Temperature, &
        k_B   => CONSTANT % BOLTZMANN, &
        m_b   => CONSTANT % ATOMIC_MASS_UNIT )

    select type ( FA => MW % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_NR ( )

    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    MeanMolecularWeight  =  k_B / ( ( Gamma - 1.0_KDR ) * C_V * m_b )

    call F % SetAdiabaticIndex ( Gamma )
    call F % SetFiducialBaryonDensity ( N_0 )
    call F % SetFiducialTemperature ( T_0 )
    call F % SetMeanMolecularWeight ( MeanMolecularWeight )

    associate &
      (   N => F % Value ( :, F % COMOVING_DENSITY ), &
        V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
        V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
          T => F % Value ( :, F % TEMPERATURE ), &
        C_S => F % Value ( :, F % SOUND_SPEED ) )

    N  =  N_0
    T  =  T_0

    V_1  =  0.0_KDR
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR

    call F % ComputeFromPrimitiveTemperature ( G )

    SoundSpeed = maxval ( C_S )

    !-- Module variable for accessibility in ApplySources_Fluid below
    Fluid => FA % Fluid_P_NR ( )

    end associate !-- N, etc.
    end select !-- PS
    end select !-- FA
    end associate !-- Gamma, etc.
    nullify ( F, G )

  end subroutine SetFluid


end module MarshakWave_Form
