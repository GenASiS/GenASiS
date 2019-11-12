module MarshakWave_Form

  !-- Vaytet et al. 2011

  use GenASiS
  use InteractionsExamples_ASC__Form
  use InteractionsExamples_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( RadiationBoxForm ) :: MarshakWaveForm
    real ( KDR ) :: &
      BoxLength, &
      AdiabaticIndex, &
      SpecificHeatCapacity, &  !-- per unit mass
      MassDensity, &
      Temperature, &
      TemperatureInner, &
      SpecificOpacity, &
      SpecificOpacityFloor, &
      EnergyMax, &
      SoundSpeed
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
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
      SetProblem

      private :: &
        InitializeInteractions, &
        SetFluid, &
        SetRadiation

    character ( LDF ), private :: &
      InteractionsType

contains


  subroutine Initialize_MW ( MW, MomentsType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    if ( MW % Type == '' ) &
      MW % Type = 'a MarshakWave'

    call InitializeRadiationBox ( MW, MomentsType, Name )
    call SetProblem ( MW, MomentsType )

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
      FinishTime, &
      MaxEnergy, &
      MinWidthEnergy, &
      EnergyScale


    !-- FinishTime

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND


    !-- Position space parameters

    allocate ( MW % BoundaryConditionsFace ( 3 ) )
    associate ( BCF => MW % BoundaryConditionsFace )
    do iD = 1, 3
      call BCF ( iD ) % Initialize ( [ 'INFLOW', 'INFLOW' ] )     
    end do
    end associate !-- BCF

    MW % BoxLength  =  25.0_KDR  *  UNIT % CENTIMETER
    call PROGRAM_HEADER % GetParameter ( MW % BoxLength, 'BoxLength' )

    MW % MinCoordinate  =  0.0_KDR
    MW % MaxCoordinate  =  MW % BoxLength


    !-- Momentum space parameters

    !-- Geometric spacing
    MinWidthEnergy  =  3.0e1_KDR  *  UNIT % KELVIN
    MaxEnergy       =  1.0e4_KDR  *  UNIT % KELVIN

    !-- Compactified spacing
    EnergyScale     =  1.0e3_KDR  *  UNIT % KELVIN


    !-- Units

    MW % Units % Time    =  UNIT % SECOND
    MW % Units % Length  =  UNIT % CENTIMETER

    MW % Units % Coordinate_PS  =  UNIT % CENTIMETER

    MW % Units % Coordinate_MS        =  UNIT % IDENTITY
    MW % Units % Coordinate_MS ( 1 )  =  UNIT % ELECTRON_VOLT

    MW % Units % BaryonMass     =  UNIT % GRAM
    MW % Units % NumberDensity  =  UNIT % MOLE  *  UNIT % CENTIMETER ** (-3)
    MW % Units % EnergyDensity  =  UNIT % ERG   *  UNIT % CENTIMETER ** (-3)
    MW % Units % Temperature    =  UNIT % KELVIN

    MW % Units % Velocity_U  =  UNIT % CENTIMETER  /  UNIT % SECOND 

    MW % Units % MomentumDensity_U &
      =  UNIT % GRAM  *  UNIT % CENTIMETER ** (-2)  /  UNIT % SECOND
    MW % Units % MomentumDensity_D &
      =  UNIT % GRAM  *  UNIT % CENTIMETER ** (-2)  /  UNIT % SECOND


    !-- Initialization

    call MW % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'PHOTONS' ], &
             MomentsType = MomentsType, &
             Name = Name, &
             MinCoordinateOption = MW % MinCoordinate, &
             MaxCoordinateOption = MW % MaxCoordinate, &
             FinishTimeOption = FinishTime, &
             MinWidthEnergyOption = MinWidthEnergy, &
             MaxEnergyOption = MaxEnergy, &
             EnergyScaleOption = EnergyScale )


  end subroutine InitializeRadiationBox


  subroutine SetProblem ( MW, MomentsType )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in )  :: &
      MomentsType

    real ( KDR ) :: &
      MeanFreePath, &
      OpticalDepth, &
      DynamicalTime, &
      DiffusionTime


    !-- Parameters

    associate &
      ( Gamma  => MW % AdiabaticIndex, &
        C_V    => MW % SpecificHeatCapacity, &
        Rho_0  => MW % MassDensity, &
        T_0    => MW % Temperature, &
        T_I    => MW % TemperatureInner, &
        Kappa  => MW % SpecificOpacity, &
        Kappa_Min => MW % SpecificOpacityFloor )

    Gamma      =  1.4_KDR
    C_V        =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    Rho_0      =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_0        =  3.0e2_KDR   *  UNIT % KELVIN
    T_I        =  1.0e3_KDR   *  UNIT % KELVIN
    Kappa      =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    Kappa_Min  =  10.0_KDR    *  UNIT % CENTIMETER ** 2 / UNIT % GRAM

    call PROGRAM_HEADER % GetParameter ( Gamma,     'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,       'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( Rho_0,     'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,       'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,       'TemperatureInner' )
    call PROGRAM_HEADER % GetParameter ( Kappa,     'SpecificOpacity' )
    call PROGRAM_HEADER % GetParameter ( Kappa_Min, 'SpecificOpacityFloor' )


    !-- Initialization

    call InitializeInteractions ( MW, MomentsType )
    call SetFluid ( MW )
    call SetRadiation ( MW )

    associate &
      ( L      => MW % BoxLength, &
        c_s    => MW % SoundSpeed, &
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
    call Show ( L, 'L' )
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
    case ( 'MARSHAK_WAVE_VAYTET_2', 'MARSHAK_WAVE_VAYTET_3' )
      call Show ( MW % EnergyMax, UNIT % ELECTRON_VOLT, 'EnergyMax' )
    end select !-- InteractionsType


    !-- Cleanup

    end associate !-- Mu, etc.
    end associate !-- Gamma, etc.

  end subroutine SetProblem


  subroutine InitializeInteractions ( MW, MomentsType )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in ) :: &
      MomentsType

    InteractionsType = 'MARSHAK_WAVE_VAYTET_1'
    call PROGRAM_HEADER % GetParameter ( InteractionsType, 'InteractionsType' )

    select case ( trim ( InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_2', 'MARSHAK_WAVE_VAYTET_3' )
      MW % EnergyMax  =  0.620_KDR * UNIT % ELECTRON_VOLT
      call PROGRAM_HEADER % GetParameter ( MW % EnergyMax, 'EnergyMax' )
    end select !-- InteractionsType

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
               MomentsType = MomentsType, Units = MW % Units )
      select case ( trim ( InteractionsType ) )
      case ( 'MARSHAK_WAVE_VAYTET_1' )
        call IA % Set_MWV_Grey &
               ( FA, SpecificOpacity = MW % SpecificOpacity )
      case ( 'MARSHAK_WAVE_VAYTET_2' )
        call IA % Set_MWV_Grey &
               ( FA, SpecificOpacity = MW % SpecificOpacity, &
                 EnergyMax = MW % EnergyMax )
      case ( 'MARSHAK_WAVE_VAYTET_3' )
        call IA % Set_MWV_Grey &
               ( FA, SpecificOpacity = MW % SpecificOpacity, &
                 EnergyMax = MW % EnergyMax, &
                 TemperatureScale = MW % Temperature )
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
             ( MS, InteractionsType = InteractionsType, Units = MW % Units )
      select case ( trim ( InteractionsType ) )
      case ( 'MARSHAK_WAVE_VAYTET_1' )
        call IB % Set_MWV_Spectral &
               ( FA, SpecificOpacity = MW % SpecificOpacity )
      case ( 'MARSHAK_WAVE_VAYTET_2' )
        call IB % Set_MWV_Spectral &
               ( FA, SpecificOpacity = MW % SpecificOpacity, &
                 SpecificOpacityFloor = MW % SpecificOpacityFloor, &
                 EnergyMax = MW % EnergyMax )
      case ( 'MARSHAK_WAVE_VAYTET_3' )
        call IB % Set_MWV_Spectral &
               ( FA, SpecificOpacity = MW % SpecificOpacity, &
                 SpecificOpacityFloor = MW % SpecificOpacityFloor, &
                 EnergyMax = MW % EnergyMax, &
                 TemperatureScale = MW % Temperature )
      end select !-- InteractionsType
      call RMB % SetInteractions ( IB )
      end select !-- IB

      end select !-- MS
      end select !-- FA
      end select !-- RMB

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
      ( Gamma => MW % AdiabaticIndex, &
        C_V   => MW % SpecificHeatCapacity, &
        Rho_0 => MW % MassDensity, &
        T_0   => MW % Temperature )

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

    MW % SoundSpeed = maxval ( C_S )

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
      ( T_0 => MW % Temperature, &
        T_I => MW % TemperatureInner, &
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

      where ( X  <  MW % MinCoordinate ( 1 ) &
              .or. Y  <  MW % MinCoordinate ( 2 ) &
              .or. Z  <  MW % MinCoordinate ( 3 ) )
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
          where ( X  <  MW % MinCoordinate ( 1 ) &
                  .or. Y  <  MW % MinCoordinate ( 2 ) &
                  .or. Z  <  MW % MinCoordinate ( 3 ) )
            J    =  J_I ( iS )
            H_1  =  0.0_KDR
            H_2  =  0.0_KDR
            H_3  =  0.0_KDR
          end where
          where ( X  >  MW % MaxCoordinate ( 1 ) &
                 .or. Y  >  MW % MaxCoordinate ( 2 ) &
                 .or. Z  >  MW % MaxCoordinate ( 3 ) )
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
