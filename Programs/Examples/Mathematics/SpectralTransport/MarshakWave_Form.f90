module MarshakWave_Form

  !-- Vaytet et al. 2011

  use Basics
  use Mathematics
  use Fluid_P_NR__Form
  use Sources_F__Form
  use Fluid_ASC__Form
  use SetPlanckSpectrum_Command
  use RadiationMoments_Form
  use PhotonMoments_S__Form
  use Sources_RM__Form
  use RadiationMoments_ASC__Form
  use Interactions_Template
  use Interactions_MWV_1_S__Form
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
    procedure, private, pass :: &
      PrepareStep_MS
    procedure, private, pass :: &
      PrepareStep_PS
  end type MarshakWaveForm

    private :: &
      SetRadiation, &
      SetFluid, &
      SetInteractions, &
      IntegrateSources, &
      ApplySources_Fluid

      private :: &
        ComputeFluidSource_G_S_Radiation_Kernel, &
        ApplySources_Fluid_Kernel

    real ( KDR ), private :: &
      AdiabaticIndex, &
      SpecificHeatCapacity, &
      MeanMolecularWeight, &
      MassDensity, &
      Temperature, &
      TemperatureInner, &
      SpecificOpacity, &
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
    class ( Interactions_BSLL_ASC_CSLD_Form ), pointer :: &
      InteractionsBundle => null ( )

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
      MeanFreePath, &
      OpticalDepth, &
      DynamicalTime, &
      DiffusionTime, &
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
      CoordinateUnit_PS, &
      CoordinateUnit_MS, &
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
        T_I    => TemperatureInner, &
        Kappa  => SpecificOpacity )

    L      =  20.0_KDR    *  UNIT % CENTIMETER
    Gamma  =  1.4_KDR
    C_V    =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    N_0    =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_0    =  3.0e2_KDR   *  UNIT % KELVIN
    T_I    =  1.0e3_KDR   *  UNIT % KELVIN
    Kappa  =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    EnergyScale  =  sqrt ( T_0 * T_I )

    call PROGRAM_HEADER % GetParameter ( L,     'BoxLength' )
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,   'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( N_0,   'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,   'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,   'TemperatureInner' )
    call PROGRAM_HEADER % GetParameter ( Kappa, 'SpecificOpacity' )
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

    CoordinateUnit_PS  =  UNIT % CENTIMETER

    MinCoordinate  =  0.0_KDR
    MaxCoordinate  =  BoxLength

    call PS % CreateChart &
           ( CoordinateUnitOption = CoordinateUnit_PS, &
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

    CoordinateUnit_MS = UNIT % IDENTITY
    CoordinateUnit_MS ( 1 ) = UNIT % ELECTRON_VOLT

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    nEnergyCells = 20
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit_MS, &
             ScaleOption = Scale, &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Prepare for Currents

    allocate ( MW % TimeStepLabel ( 2 ) )
    MW % TimeStepLabel ( 1 )  =  'Radiation'
    MW % TimeStepLabel ( 2 )  =  'Fluid'

    TimeUnit = UNIT % SECOND

    VelocityUnit ( 1 )     =  CoordinateUnit_PS ( 1 ) / TimeUnit 
    VelocityUnit ( 2 )     =  CoordinateUnit_PS ( 2 ) / TimeUnit
    VelocityUnit ( 3 )     =  CoordinateUnit_PS ( 3 ) / TimeUnit
    MassDensityUnit        =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit      =  UNIT % MASS_DENSITY_CGS  &
                              *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumDensity_U_Unit =  EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
    MomentumDensity_D_Unit =  EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
    TemperatureUnit        =  UNIT % KELVIN

    MassUnit  =  MassDensityUnit  *  CoordinateUnit_PS ( 1 )
    if ( PS % nDimensions > 1 ) &
      MassUnit  =  MassUnit  *  CoordinateUnit_PS ( 2 )
    if ( PS % nDimensions > 2 ) &
      MassUnit  =  MassUnit  *  CoordinateUnit_PS ( 3 )

    EnergyUnit           =  MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumUnit         =  MassUnit  *  UNIT % SPEED_OF_LIGHT
    AngularMomentumUnit  =  CoordinateUnit_PS ( 1 ) &
                            *  MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    !-- Radiation

    allocate ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
               MW % Current_BSLL_ASC_CSLD )
    select type ( RMB => MW % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )
    call RMB % Initialize &
           ( MS, 'PHOTONS_SPECTRAL', &
             MomentumDensity_U_UnitOption &
               = MomentumDensity_U_Unit / CoordinateUnit_MS ( 1 ) ** 3, &
             MomentumDensity_D_UnitOption &
               = MomentumDensity_D_Unit / CoordinateUnit_MS ( 1 ) ** 3, &
             EnergyDensityUnitOption &
               = EnergyDensityUnit / CoordinateUnit_MS ( 1 ) ** 3, &
             TimeUnitOption = TimeUnit )

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
             AngularMomentumUnitOption = AngularMomentumUnit, &
             TimeUnitOption = TimeUnit )

    !-- Interactions

    allocate ( MW % Interactions_BSLL_ASC_CSLD )
    associate ( IB => MW % Interactions_BSLL_ASC_CSLD )
    call IB % Initialize &
           ( MS, InteractionsType, &
             LengthUnitOption = CoordinateUnit_PS ( 1 ), &
             EnergyDensityUnitOption = EnergyDensityUnit )
    call RMB % SetInteractions ( IB )
    end associate !-- IB

    !-- Step

    allocate ( Step_RK2_C_BSLL_ASC_CSLD_Form :: MW % Step_MS )
    select type ( S_MS => MW % Step_MS )
    class is ( Step_RK2_C_BSLL_ASC_CSLD_Form )
    call S_MS % Initialize ( RMB, Name )
    S_MS % ApplyRelaxation_F % Pointer => ApplyRelaxation_RM
    end select !-- S_MS

    allocate ( Step_RK2_C_ASC_Form :: MW % Step_PS )
    select type ( S_PS => MW % Step_PS )
    class is ( Step_RK2_C_ASC_Form )
    call S_PS % Initialize ( FA, Name )
    S_PS % ApplySources % Pointer  =>  ApplySources_Fluid
    end select !-- S_PS

    !-- Initial conditions

    call SetRadiation ( MW )
    call SetFluid ( MW )
    call SetInteractions ( MW )

    associate &
      ( Mu     => MeanMolecularWeight, & 
        c_s    => SoundSpeed, &
        t_Dyn  => DynamicalTime, &
        Lambda => MeanFreePath, &
        Tau    => OpticalDepth, &
        t_Diff => DiffusionTime, &
        c      => CONSTANT % SPEED_OF_LIGHT &
      )

    t_Dyn   =  L / c_s

    Lambda  =  1.0_KDR / ( N_0 * Kappa )
    Tau     =  L / Lambda
    t_Diff  =  L * Tau / c

    call Show ( 'MarshakWave parameters' )
    call Show ( Gamma, 'Gamma' )
    call Show ( C_V, UNIT % ERG / UNIT % KELVIN / UNIT % GRAM, 'C_V' )
    call Show ( Mu, 'Mu' )
    call Show ( N_0, UNIT % MASS_DENSITY_CGS, 'N_0' )
    call Show ( T_0, UNIT % KELVIN, 'T_0' )
    call Show ( c_s, UNIT % SPEED_OF_LIGHT, 'c_s' )
    call Show ( t_Dyn, UNIT % SECOND, 't_Dyn' )
    call Show ( T_I, UNIT % KELVIN, 'T_I' )
    call Show ( Kappa, UNIT % CENTIMETER ** 2 / UNIT % GRAM, 'Kappa' )
    call Show ( Lambda, UNIT % CENTIMETER, 'Lambda' )
    call Show ( Tau, UNIT % IDENTITY, 'Tau' )
    call Show ( t_Diff, UNIT % SECOND, 't_Diff' )
    call Show ( InteractionsType, 'InteractionsType' )

    select case ( trim ( InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_2_GREY', 'MARSHAK_WAVE_VAYTET_3_GREY' )
      call Show ( EnergyMax, UNIT % ELECTRON_VOLT, 'EnergyMax' )
    end select !-- InteractionsType

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


  subroutine PrepareStep_MS ( I )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( PhotonMoments_S_Form ), pointer :: &
      Radiation
    class ( InteractionsTemplate ), pointer :: &
      Interactions

!call Show ( '>>> Prepare_MS' )
    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    do iF = 1, MS % nFibers
      Radiation    =>  RadiationBundle % PhotonMoments_S ( iF )
      Interactions =>  InteractionsBundle % Interactions ( iF )
      call Clear ( Radiation % Sources % Value )
!if ( iF == 1 ) then
!call Show ( iF, '>>> iF' )
!call Show ( Radiation % Value ( :, Radiation % COMOVING_ENERGY_EQ ), '>>> J_Eq pre' )
!end if
      call Interactions % Compute ( Radiation )
!if ( iF == 1 ) then
!call Show ( Radiation % Value ( :, Radiation % COMOVING_ENERGY_EQ ), '>>> J_Eq' )
!call Show ( Radiation % Value ( :, Radiation % COMOVING_ENERGY ), '>>> J' )
!end if
    end do !-- iF

    end select !-- MS

    nullify ( Radiation, Interactions )

  end subroutine PrepareStep_MS


  subroutine PrepareStep_PS ( I )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      I

    class ( RadiationMomentsForm ), pointer :: &
      RMEI

!call Show ( '>>> Prepare_PS' )

    call IntegrateSources ( I )

    RMEI => RadiationBundle % EnergyIntegral % RadiationMoments ( )

    select type ( SR => RMEI % Sources )
    class is ( Sources_RM_Form )

    select type ( SF => Fluid % Sources )
    class is ( Sources_F_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    call Clear ( SF % Value )

    call ComputeFluidSource_G_S_Radiation_Kernel &
           ( SF % Value ( :, SF % RADIATION_G ), & 
             SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
             CSL % IsProperCell, &
             SR % Value ( :, SR % INTERACTIONS_J ), &
             SR % Value ( :, SR % INTERACTIONS_H_D ( 1 ) ) )

    end select !-- CSL
    end select !-- PS
    end select !-- SF
    end select !-- SR

    nullify ( RMEI )

  end subroutine PrepareStep_PS


  subroutine SetRadiation ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iS     !-- iSection
    real ( KDR ), dimension ( : ), allocatable :: &
      J_0, &
      J_I
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( RadiationMomentsForm ), pointer :: &
      RS
    class ( PhotonMoments_S_Form ), pointer :: &
      RM

    select type ( MS => MW % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    G => MS % Base_CSLD % Geometry ( )

    select type ( RMB => MW % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    associate &
      ( T_0  =>  Temperature, &
        T_I  =>  TemperatureInner, &
        E    =>  RMB % Energy, &
        X    =>  G % Value ( :, G % CENTER ( 1 ) ), &
        Y    =>  G % Value ( :, G % CENTER ( 2 ) ), &
        Z    =>  G % Value ( :, G % CENTER ( 3 ) ) )

    !-- Proper cells

    do iF = 1, MS % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )
      RM => RMB % PhotonMoments_S ( iF )
      associate &
        ( J     =>  RM % Value ( :, RM % COMOVING_ENERGY ), &
          J_Eq  =>  RM % Value ( :, RM % COMOVING_ENERGY_EQ ), &
          H_1   =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
          H_2   =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
          H_3   =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ) )

      call SetPlanckSpectrum ( E, T_0, J )
      call SetPlanckSpectrum ( E, T_0, J_Eq )

      H_1  =  0.0_KDR
      H_2  =  0.0_KDR
      H_3  =  0.0_KDR

      call RM % ComputeFromPrimitive ( iBC, G )

      end associate !-- J, etc.
      end associate !-- iBC
    end do !-- iF

    call RMB % LoadSections ( )

    !-- Boundary cells

    allocate &
      ( J_0 ( MS % nSections ), &
        J_I ( MS % nSections ) )

    call SetPlanckSpectrum ( E, T_I, J_I )
    call SetPlanckSpectrum ( E, T_0, J_0 )

    do iS = 1, MS % nSections
      select type ( RSA => RMB % Section % Atlas ( iS ) % Element )
      class is ( RadiationMoments_ASC_Form )
        RS => RSA % RadiationMoments ( )
        associate &
          ( J     =>  RS % Value ( :, RS % COMOVING_ENERGY ), &
            J_Eq  =>  RS % Value ( :, RS % COMOVING_ENERGY_EQ ), &
            H_1   =>  RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 1 ) ), &
            H_2   =>  RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 2 ) ), &
            H_3   =>  RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 3 ) ) )
        where ( X < MinCoordinate ( 1 ) .or. Y < MinCoordinate ( 2 ) &
               .or. Z < MinCoordinate ( 3 ) )
          J     =  J_I ( iS )
          J_Eq  =  J_I ( iS )
          H_1   =  0.0_KDR
          H_2   =  0.0_KDR
          H_3   =  0.0_KDR
        end where
        where ( X > MinCoordinate ( 1 ) .or. Y > MinCoordinate ( 2 ) &
               .or. Z > MinCoordinate ( 3 ) )
          J     =  J_0 ( iS )
          J_Eq  =  J_0 ( iS )
          H_1   =  0.0_KDR
          H_2   =  0.0_KDR
          H_3   =  0.0_KDR
        end where
        end associate !-- J, etc.
        call RS % ComputeFromPrimitive ( G )
      end select
    end do

    deallocate ( J_0, J_I )
    
    !-- Module variable for accessibility
    RadiationBundle => RMB
  
    end associate !-- T_0, etc.
    end select !-- RMB
    end select !-- MS

    nullify ( G, RS, RM )

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


  subroutine SetInteractions ( MW )

    type ( MarshakWaveForm ), intent ( inout ), target :: &
      MW

    integer ( KDI ) :: &
      iF  !-- iFiber
    class ( PhotonMoments_S_Form ), pointer :: &
      R
    class ( InteractionsTemplate ), pointer :: &
      I

    associate &
      ( IB   =>  MW % Interactions_BSLL_ASC_CSLD, &
        RMB  =>  RadiationBundle, &
        F    =>  Fluid )

    select type ( MS => MW % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    do iF = 1, MS % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )
      I  =>   IB % Interactions ( iF )
      R  =>  RMB % PhotonMoments_S ( iF )
      select type ( I )
      type is ( Interactions_MWV_1_S_Form )
        call I % Set ( R, F, RMB % Energy, SpecificOpacity, iBC )
    ! type is ( Interactions_MWV_2_G_Form )
    !   call I % Set ( R, F, SpecificOpacity, EnergyMax )
    ! type is ( Interactions_MWV_3_G_Form )
    !   call I % Set ( R, F, SpecificOpacity, EnergyMax, Temperature )
    class default
      call Show ( 'Interactions type not recognized', CONSOLE % ERROR )
      call Show ( 'MarshakWave_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetInteractions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I
      end associate !-- iBC
    end do !-- iF

    !-- Module variable for accessibility in PrepareStep_MS
    InteractionsBundle => IB

    end select !-- MS
    end associate !-- IB, etc.
    nullify ( R, I )

  end subroutine SetInteractions


  subroutine IntegrateSources ( MW )

    type ( MarshakWaveForm ), intent ( inout ), target :: &
      MW

    integer ( KDI ) :: &
      iI, &  !-- iIntegral
      iF, &  !-- iFiber
      nIntegrals
    real ( KDR ), dimension ( : ), allocatable :: &
      Integral
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      Integrand
    type ( VolumeIntegralForm ) :: &
      VI
    class ( RadiationMomentsForm ), pointer :: &
      RMEI, &
      RMF

    associate &
      ( IB => InteractionsBundle, &
        RMB => RadiationBundle, &
        MS => RadiationBundle % Bundle_SLL_ASC_CSLD )

    RMEI => RMB % EnergyIntegral % RadiationMoments ( )

    nIntegrals = 2

    allocate ( Integral  ( nIntegrals ) )
    allocate ( Integrand ( nIntegrals ) )
    do iI = 1, nIntegrals
      call Integrand ( iI ) % Initialize ( RMB % nEnergyValues )
    end do !-- iC

    do iF = 1, MS % nFibers
      associate &
        ( iBC => MS % iaBaseCell ( iF ), &
          CF  => MS % Fiber_CSLL )

      RMF => RMB % RadiationMoments ( iF )
      select type ( SF => RMF % Sources )
      class is ( Sources_RM_Form )
        Integrand ( 1 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_J )
        Integrand ( 2 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_H_D ( 1 ) )
      end select !-- SF

      call VI % Compute ( CF, Integrand, Integral )

      select type ( SEI => RMEI % Sources )
      class is ( Sources_RM_Form )
        SEI % Value ( iBC, SEI % INTERACTIONS_J ) = Integral ( 1 ) 
        SEI % Value ( iBC, SEI % INTERACTIONS_H_D ) = Integral ( 1 ) 
      end select !-- SEI
      
      end associate !-- iBC, etc.
    end do !-- iF

    end associate !-- IB, etc.

    nullify ( RMEI, RMF )
 
  end subroutine IntegrateSources


  subroutine ApplySources_Fluid &
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

    integer ( KDI ) :: &
      iEnergy, &
      iMomentum_1

    select type ( SF => Sources_F )
    class is ( Sources_F_Form )

    select type ( F => Fluid )
    class is ( Fluid_P_NR_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )

    call ApplySources_Fluid_Kernel &
           ( Increment % Value ( :, iEnergy ), &
             Increment % Value ( :, iMomentum_1 ), &
             Chart % IsProperCell, &
             SF % Value ( :, SF % RADIATION_G ), &
             SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
             TimeStep )

    end select !-- Chart
    end select !-- F
    end select !-- SF

  end subroutine ApplySources_Fluid

  
  subroutine ComputeFluidSource_G_S_Radiation_Kernel &
               ( FS_R_G, FS_R_S_1, IsProperCell, Q, A_1 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_G, &
      FS_R_S_1
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Q, &
      A_1

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( FS_R_G )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      FS_R_G ( iV )    =    FS_R_G ( iV )  -    Q ( iV ) 
      FS_R_S_1 ( iV )  =  FS_R_S_1 ( iV )  -  A_1 ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_G_S_Radiation_Kernel


  subroutine ApplySources_Fluid_Kernel &
               ( K_G, K_S_1, IsProperCell, FS_R_G, FS_R_S_1, dt )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_G, &
      K_S_1
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FS_R_G, &
      FS_R_S_1
    real ( KDR ), intent ( in ) :: &
      dt

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( K_G )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      K_G ( iV )    =  K_G ( iV )    +  FS_R_G ( iV )    *  dt
      K_S_1 ( iV )  =  K_S_1 ( iV )  +  FS_R_S_1 ( iV )  *  dt
    end do
    !$OMP end parallel do

  end subroutine ApplySources_Fluid_Kernel


end module MarshakWave_Form
