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
      SetFluid, &
      SetRadiation

    real ( KDR ), private :: &
      BoxLength, &
      AdiabaticIndex, &
      SpecificHeatCapacity, &  !-- per unit mass
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

    real ( KDR ) :: &
      MeanFreePath, &
      OpticalDepth, &
      DynamicalTime, &
      DiffusionTime

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
    call SetRadiation ( MW )

    associate &
      ( c_s    => SoundSpeed, &
        t_Dyn  => DynamicalTime, &
        Lambda => MeanFreePath, &
        Tau    => OpticalDepth, &
        t_Diff => DiffusionTime, &
        c      => CONSTANT % SPEED_OF_LIGHT )

    t_Dyn   =  L / c_s

    Lambda  =  1.0_KDR / ( Rho_0 * Kappa )
    Tau     =  L / Lambda
    t_Diff  =  L * Tau / c

    call Show ( 'MarshakWave parameters' )
    call Show ( Gamma, 'Gamma' )
    call Show ( C_V, UNIT % ERG / UNIT % KELVIN / UNIT % GRAM, 'C_V' )
    call Show ( Rho_0, UNIT % MASS_DENSITY_CGS, 'Rho_0' )
    call Show ( T_0, UNIT % KELVIN, 'T_0' )
    call Show ( c_s, UNIT % CENTIMETER / UNIT % SECOND, 'c_s' )
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


    !-- Cleanup

    end associate !-- Mu, etc.
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

    type ( MeasuredValueForm ) :: &
      LengthUnit, &
      EnergyDensityUnit

    LengthUnit         =  UNIT % CENTIMETER
    EnergyDensityUnit  =  UNIT % ERG  *  UNIT % CENTIMETER ** (-3)

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
               MomentsType = MomentsType, LengthUnitOption = LengthUnit, &
               EnergyDensityUnitOption = EnergyDensityUnit )
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
      call IB % Initialize &
             ( MS, InteractionsType = InteractionsType, &
               LengthUnitOption = LengthUnit, &
               EnergyDensityUnitOption = EnergyDensityUnit )
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

    E_0  =  C_V * Rho_0 * T_0
    P_0  =  ( Gamma - 1.0_KDR ) * E_0

    m_b  =  CONSTANT % ATOMIC_MASS_UNIT
    N_0  =  Rho_0 / m_b

    call F % SetAdiabaticIndex ( Gamma )
    call F % SetSpecificHeatVolume ( C_V * m_b )
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


  subroutine SetRadiation ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iS     !-- iSection
    real ( KDR ), dimension ( : ), allocatable :: &
      J_0, &
      J_I
    real ( KDR ) :: &
      a
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RS

    !-- Fluid and Geometry

    select type ( I => MW % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )
      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )
        F => FA % Fluid_P_I ( )
      end select !-- FA
      select type ( PS => I % PositionSpace )
      class is ( Atlas_SC_Form )
        G => PS % Geometry ( )
      end select !-- PS
    end select !-- I

    !-- Radiation

    associate &
      ( T_0 => Temperature, &
        T_I => TemperatureInner, &
          X => G % Value ( :, G % CENTER_U ( 1 ) ), &
          Y => G % Value ( :, G % CENTER_U ( 2 ) ), &
          Z => G % Value ( :, G % CENTER_U ( 3 ) ) )

    select type ( I => MW % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )
      RM => RMA % RadiationMoments ( )

      associate &
        (   J => RM % Value ( :, RM % COMOVING_ENERGY ), &
          H_1 => RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
          H_2 => RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
          H_3 => RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ) )
      
      a  =  4.0_KDR * CONSTANT % STEFAN_BOLTZMANN

      J  =  a  *  T_0 ** 4

      H_1  =  0.0_KDR
      H_2  =  0.0_KDR
      H_3  =  0.0_KDR

      where ( X < MinCoordinate ( 1 ) .or. Y < MinCoordinate ( 2 ) &
              .or. Z < MinCoordinate ( 3 ) )
        J  =  a  *  T_I ** 4
      end where

      call RM % ComputeFromPrimitive ( G )

      ! !-- Module variable for accessibility
      ! Radiation => RA % PhotonMoments_G ( )

      end associate !-- J, etc.
      end select !-- RMA

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )

      associate ( E => RMB % Energy )

      !-- Proper cells

      do iF = 1, MS % nFibers
        associate ( iBC => MS % iaBaseCell ( iF ) )
        RM => RMB % RadiationMoments ( iF )
        associate &
          ( J    =>  RM % Value ( :, RM % COMOVING_ENERGY ), &
            H_1  =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
            H_2  =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
            H_3  =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ) )

        call SetPlanckSpectrum ( E, T_0, J )

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
              H_1   =>  RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 1 ) ), &
              H_2   =>  RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 2 ) ), &
              H_3   =>  RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 3 ) ) )
          where ( X < MinCoordinate ( 1 ) .or. Y < MinCoordinate ( 2 ) &
                 .or. Z < MinCoordinate ( 3 ) )
            J    =  J_I ( iS )
            H_1  =  0.0_KDR
            H_2  =  0.0_KDR
            H_3  =  0.0_KDR
          end where
          where ( X > MaxCoordinate ( 1 ) .or. Y > MaxCoordinate ( 2 ) &
                 .or. Z > MaxCoordinate ( 3 ) )
            J    =  J_0 ( iS )
            H_1  =  0.0_KDR
            H_2  =  0.0_KDR
            H_3  =  0.0_KDR
          end where
          end associate !-- J, etc.
          call RS % ComputeFromPrimitive ( G )
        end select
      end do !-- iS

      deallocate ( J_0, J_I )
    
!     !-- Module variable for accessibility
!     RadiationBundle => RMB
  
      end associate !-- E, etc.
      end select !-- MS
      end select !-- RMB

    end select !-- I
    end associate !-- T_0, etc.
    nullify ( G, F, RM, RS )

  end subroutine SetRadiation


end module MarshakWave_Form
