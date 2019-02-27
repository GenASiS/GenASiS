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
      InitializeInteractions, &
      SetFluid

    real ( KDR ), private :: &
      BoxLength, &
      AdiabaticIndex, &
      SpecificHeatCapacity, &
      MassDensity, &
      Temperature, &
      TemperatureInner, &
      SpecificOpacity, &
      SpecificOpacityFloor, &
      EnergyMax, &
      SoundSpeed, &
      FinishTime
    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDF ) :: &
      InteractionsType

!    class ( MarshakWaveForm ), private, pointer :: &
!      MarshakWave => null ( )

contains


  subroutine Initialize_MW ( MW, MomentsType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    if ( MW % Type == '' ) &
      MW % Type = 'a MarshakWave'

!    MarshakWave => MW


    !-- Parameters

    associate &
      ( L      => BoxLength, &
        Gamma  => AdiabaticIndex, &
        C_V    => SpecificHeatCapacity, &
        Rho_0  => MassDensity, &
        T_0    => Temperature, &
        T_I    => TemperatureInner, &
        Kappa  => SpecificOpacity, &
        Kappa_Min => SpecificOpacityFloor )

    L          =  20.0_KDR    *  UNIT % CENTIMETER
    Gamma      =  1.4_KDR
    C_V        =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    Rho_0      =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_0        =  3.0e2_KDR   *  UNIT % KELVIN
    T_I        =  1.0e3_KDR   *  UNIT % KELVIN
    Kappa      =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    Kappa_Min  =  10.0_KDR    *  UNIT % CENTIMETER ** 2 / UNIT % GRAM

    call PROGRAM_HEADER % GetParameter ( L,         'BoxLength' )
    call PROGRAM_HEADER % GetParameter ( Gamma,     'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,       'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( Rho_0,     'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,       'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,       'TemperatureInner' )
    call PROGRAM_HEADER % GetParameter ( Kappa,     'SpecificOpacity' )
    call PROGRAM_HEADER % GetParameter ( Kappa_Min, 'SpecificOpacityFloor' )

    InteractionsType = 'MARSHAK_WAVE_VAYTET_1'
    call PROGRAM_HEADER % GetParameter ( InteractionsType, 'InteractionsType' )

    select case ( trim ( InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_2', 'MARSHAK_WAVE_VAYTET_3' )
      EnergyMax  =  0.620_KDR * UNIT % ELECTRON_VOLT
      call PROGRAM_HEADER % GetParameter ( EnergyMax, 'EnergyMax' )
    end select !-- InteractionsType

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND


    !-- Initialization

    call InitializeRadiationBox ( MW, MomentsType, Name )
    call InitializeInteractions ( MW, MomentsType )
    call SetFluid ( MW )


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
    call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )


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
             RadiationType = [ 'PHOTONS' ], &
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


  subroutine InitializeInteractions ( MW, MomentsType )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in ) :: &
      MomentsType

    select type ( I => MW % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )

      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )

      select type ( PS => I % PositionSpace )
      class is ( Atlas_SC_Form )

      allocate ( InteractionsExamples_ASC_Form :: MW % Interactions_ASC )
      select type ( IA => MW % Interactions_ASC )
      class is ( InteractionsExamples_ASC_Form )
      call IA % Initialize &
             ( PS, InteractionsType = InteractionsType, &
               MomentsType = MomentsType )
      select case ( trim ( InteractionsType ) )
      case ( 'MARSHAK_WAVE_VAYTET_1' )
        call IA % Set_MWV_Grey ( FA, SpecificOpacity = SpecificOpacity )
      end select !-- InteractionsType
      call RMA % SetInteractions ( IA )
      end select !-- IA

      end select !-- PS
      end select !-- FA
      end select !-- RMA

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      class is ( RadiationMoments_BSLL_ASC_CSLD_Form )

      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )

      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )

      allocate ( InteractionsExamples_BSLL_ASC_CSLD_Form :: &
                   MW % Interactions_BSLL_ASC_CSLD )
      select type ( IB => MW % Interactions_BSLL_ASC_CSLD )
      class is ( InteractionsExamples_BSLL_ASC_CSLD_Form )
      call IB % Initialize ( MS, InteractionsType = InteractionsType )
      select case ( trim ( InteractionsType ) )
      case ( 'MARSHAK_WAVE_VAYTET_1' )
        call IB % Set_MWV_Spectral ( FA, SpecificOpacity = SpecificOpacity )
      end select !-- InteractionsType
      call RMB % SetInteractions ( IB )
      end select !-- IB

      end select !-- PS
      end select !-- FA
      end select !-- RMA

    end select !-- Integrator

  end subroutine InitializeInteractions


  subroutine SetFluid ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    real ( KDR ) :: &
      m_b, &
      N_0, &
      E_0, &
      P_0

    associate &
      ( Gamma => AdiabaticIndex, &
        C_V   => SpecificHeatCapacity, &
        Rho_0 => MassDensity, &
        T_0   => Temperature )

    select type ( I => MW % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    m_b  =  CONSTANT % ATOMIC_MASS_UNIT
    N_0  =  Rho_0 / m_b
    E_0  =  C_V * N_0 * T_0
    P_0  =  ( Gamma - 1.0_KDR ) * E_0

    call F % SetAdiabaticIndex ( Gamma )
    call F % SetSpecificHeatVolume ( C_V )
    call F % SetFiducialParameters ( N_0, P_0 )

    associate &
      (   N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
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

    call F % ComputeFromTemperature ( F % Value, G, G % Value )

    SoundSpeed = maxval ( C_S )

    ! !-- Module variable for accessibility in ApplySources_Fluid below
    ! Fluid => FA % Fluid_P_NR ( )

    end associate !-- N, etc.
    end select !-- PS
    end select !-- FA
    end select !-- I
    end associate !-- Gamma, etc.
    nullify ( F, G )

  end subroutine SetFluid


end module MarshakWave_Form
