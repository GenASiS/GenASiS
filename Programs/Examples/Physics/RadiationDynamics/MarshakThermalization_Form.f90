module MarshakThermalization_Form

  !-- A thermalization problem inspired by setup of Vaytet et al. 2011

  use GenASiS
  use Interactions_MWV_1__Form
  use Interactions_MWV_2__Form
  use Interactions_MWV_3__Form

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: MarshakThermalizationForm
    real ( KDR ) :: &  !-- Primary parameters
      BoxLength, &
      AdiabaticIndex, &
      SpecificHeatCapacity, &  !-- per unit mass
      MassDensity, &
      TemperatureFluid, &
      TemperatureRadiation, &
      SpecificOpacity, &
      SpecificOpacityMin, &
      EnergyMax
   real ( KDR ) :: &  !-- Derived parameters
      SoundSpeed, &
      DynamicalTime, &
      MeanFreePath, &
      OpticalDepth, &
      DiffusionTime
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDL ) :: &
      InteractionsType = ''
  contains
    procedure, public, pass :: &
      Initialize_MT
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type MarshakThermalizationForm

    private :: &
      InitializeUniverse, &
      SetInitial

     private :: &
       SetFluid, &
       SetRadiation

        private :: &
          SetFluidKernel, &
          SetRadiationKernel


contains


  subroutine Initialize_MT ( MT, FormalismType, Name )

    class ( MarshakThermalizationForm ), intent ( inout ), target :: &
      MT
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( MT % Type  ==  '' ) &
      MT % Type  =  'a MarshakThermalization'

    call InitializeUniverse ( MT, FormalismType, Name )
 
  end subroutine Initialize_MT


  impure elemental subroutine Finalize ( MT )

    type ( MarshakThermalizationForm ), intent ( inout ), target :: &
      MT

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( MarshakThermalizationForm ), intent ( in ) :: &
      U

    call U % Universe_R_B_Form % ShowParameters ( )

    call Show ( 'Primary parameters' )
    call Show ( U % BoxLength, &
                UNIT % CENTIMETER, &
                'BoxLength' )
    call Show ( U % AdiabaticIndex, &
                'AdiabaticIndex' )
    call Show ( U % SpecificHeatCapacity, &
                UNIT % ERG / UNIT % KELVIN / UNIT % GRAM, &
                'SpecificHeatCapacity' )
    call Show ( U % MassDensity, &
                UNIT % MASS_DENSITY_CGS, &
                'MassDensity' )
    call Show ( U % TemperatureFluid, &
                UNIT % KELVIN, &
                'TemperatureFluid' )
    call Show ( U % TemperatureRadiation, &
                UNIT % KELVIN, &
                'TemperatureRadiation' )
    call Show ( U % SpecificOpacity, &
                UNIT % CENTIMETER ** 2 / UNIT % GRAM, &
                'SpecificOpacity' )
    call Show ( U % SpecificOpacityMin, &
                UNIT % CENTIMETER ** 2 / UNIT % GRAM, &
                'SpecificOpacityMin' )
    call Show ( U % EnergyMax, &
                UNIT % ELECTRON_VOLT, &
                'EnergyMax' )
    call Show ( U % InteractionsType, &
                'InteractionsType' )
                
    call Show ( 'Derived Parameters' )
    call Show ( U % SoundSpeed, &
                UNIT % SPEED_CGS, &
                'SoundSpeed' )
    call Show ( U % DynamicalTime, &
                UNIT % SECOND, &
                'DynamicalTime' )

  end subroutine ShowParameters


  subroutine InitializeUniverse ( MT, FormalismType, Name )

    class ( MarshakThermalizationForm ), intent ( inout ), target :: &
      MT
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    integer ( KDI ) :: &
      iD

    !-- Position space parameters

    MT % BoxLength  =  25.0_KDR  *  UNIT % CENTIMETER
    call PROGRAM_HEADER % GetParameter ( MT % BoxLength, 'BoxLength' )

    MT % MinCoordinate  =  0.0_KDR
    MT % MaxCoordinate  =  MT % BoxLength

    ! !-- Momentum space parameters

    ! !-- Geometric spacing
    ! MinWidthEnergy  =  3.0e1_KDR  *  UNIT % KELVIN
    ! MaxEnergy       =  1.0e4_KDR  *  UNIT % KELVIN

    ! !-- Compactified spacing
    ! EnergyScale     =  1.0e3_KDR  *  UNIT % KELVIN


    !-- Interactions

    MT % InteractionsType = 'MARSHAK_WAVE_VAYTET_1'
    call PROGRAM_HEADER % GetParameter &
           ( MT % InteractionsType, 'InteractionsType' )

    select case ( trim ( MT % InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_1' )
      allocate ( Interactions_MWV_1_Form :: MT % Interactions_BM )
    case ( 'MARSHAK_WAVE_VAYTET_2' )
      allocate ( Interactions_MWV_2_Form :: MT % Interactions_BM )
    case ( 'MARSHAK_WAVE_VAYTET_3' )
      allocate ( Interactions_MWV_3_Form :: MT % Interactions_BM )
    case default
      call Show ( 'InteractionsType not recognized', CONSOLE % ERROR )
      call Show ( 'MarshakWave_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeUniverse', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    !-- Initialization

    call MT % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'PHOTONS' ], &
             FormalismType = FormalismType, &
             Name = Name, &
             UnitsTypeOption = 'CGS', &
             MinCoordinateOption = MT % MinCoordinate, &
             MaxCoordinateOption = MT % MaxCoordinate, &
             nCellsPositionOption = [ 128, 128, 128 ] )

    !-- Boundary conditions

    select type ( I  =>  MT % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    select type ( I  =>  MT % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    associate &
      ( R  =>  I % CurrentSet_X_1D )
    do iD  =  1, 3
      call R % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- R
    end select !-- I

    !-- Pointers

    MT % Integrator % SetInitial    =>  SetInitial
    MT % Integrator % System        =>  MT
    
  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( MT  =>  I % System )
      class is ( MarshakThermalizationForm )
    select type ( I )
      class is ( Integrator_CS_1D_CS_Form )

    !-- Parameters

    associate &
      ( Gamma      =>  MT % AdiabaticIndex, &
        C_V        =>  MT % SpecificHeatCapacity, &
        Rho_0      =>  MT % MassDensity, &
        T_F        =>  MT % TemperatureFluid, &
        T_R        =>  MT % TemperatureRadiation, &
        Kappa      =>  MT % SpecificOpacity, &
        Kappa_Min  =>  MT % SpecificOpacityMin, &
        E_Max      =>  MT % EnergyMax )

    Gamma      =  1.4_KDR
    C_V        =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    Rho_0      =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_F        =  3.0e2_KDR   *  UNIT % KELVIN
    T_R        =  1.0e3_KDR   *  UNIT % KELVIN
    Kappa      =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
!-- More diffusive
!    Kappa      =  1.0e4_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    Kappa_Min  =  10.0_KDR    *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    E_Max      =  0.620_KDR   *  UNIT % ELECTRON_VOLT

    call PROGRAM_HEADER % GetParameter ( Gamma,     'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,       'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( Rho_0,     'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_F,       'TemperatureFluid' )
    call PROGRAM_HEADER % GetParameter ( T_R,       'TemperatureRadiation' )
    call PROGRAM_HEADER % GetParameter ( Kappa,     'SpecificOpacity' )
    call PROGRAM_HEADER % GetParameter ( Kappa_Min, 'SpecificOpacityMin' )
    call PROGRAM_HEADER % GetParameter ( E_Max,     'EnergyMax' )

    !-- FinishTime

    associate ( c  =>  CONSTANT % SPEED_OF_LIGHT )
    I % T_Finish  =  100.0_KDR  *  1.0  /  ( c * Kappa * Rho_0 ) 
    end associate !-- c

    !-- Fluid

    select type ( F  =>  I % CurrentSet_X )
    class is ( Fluid_P_I_Form )
      call SetFluid ( MT, F )
      call F % SetUseInitialTemperature ( .true. )
    end select !-- F

    !-- Radiation

    select type ( I )
    class is ( Integrator_CS_1D_BM_CS_Form )

      select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )

      call SetRadiation ( MT, R )

      end select !-- R

    class default
      call Show ( 'Integrator type not recognized', CONSOLE % ERROR )
      call Show ( 'MarshakThermalizationForm', 'module', CONSOLE % ERROR )
      call Show ( 'SetInitial', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I    

    !-- Interactions

    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( PhotonMoments_G_Form )

    select type ( Intrctns  =>  MT % Interactions_BM )
    class is ( Interactions_MWV_3_Form )
       call Intrctns % SetSpecificOpacity ( MT % SpecificOpacity )
       call Intrctns % SetEnergyMax ( MT % EnergyMax )
       call Intrctns % SetTemperatureScale ( MT % TemperatureFluid )
    class is ( Interactions_MWV_2_Form )
       call Intrctns % SetSpecificOpacity ( MT % SpecificOpacity )
       call Intrctns % SetEnergyMax ( MT % EnergyMax )
    class is ( Interactions_MWV_1_Form )
       call Intrctns % SetSpecificOpacity ( MT % SpecificOpacity )
    end select !-- Intrctns
    
    end select !-- R
    end select !-- I

    !-- Cleanup

    end associate !-- Gamma, etc.
    end select !-- I
    end select !-- MT

  end subroutine SetInitial


  subroutine SetFluid ( MT, F )

    class ( MarshakThermalizationForm ), intent ( inout ) :: &
      MT
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    real ( KDR ) :: &
      N_0, &
      E_0, &
      P_0

    associate &
      (     L      =>  MT % BoxLength, &
        Gamma      =>  MT % AdiabaticIndex, &
            C_V    =>  MT % SpecificHeatCapacity, &
          Rho_0    =>  MT % MassDensity, &
            T_0    =>  MT % TemperatureFluid, &
            c_s    =>  MT % SoundSpeed, &
            t_Dyn  =>  MT % DynamicalTime, &
            m_b    =>   F % BaryonMass )
    associate &
      ( FV  =>  F % Storage_GS % Value )

    N_0  =  Rho_0 / m_b
    E_0  =  C_V * Rho_0 * T_0
    P_0  =  ( Gamma - 1.0_KDR ) * E_0

    c_s    =  sqrt ( Gamma * P_0 / Rho_0 )
    t_Dyn  =  L / c_s

    call F % SetAdiabaticIndex ( Gamma )
    call F % SetSpecificHeatVolume ( C_V * m_b )
    call F % SetFiducialParameters ( N_0, P_0 )

    call SetFluidKernel &
           ( N    =  FV ( :, F % BARYON_DENSITY_C ), &
             V_1  =  FV ( :, F % VELOCITY_U_1 ), &
             V_2  =  FV ( :, F % VELOCITY_U_2 ), &
             V_3  =  FV ( :, F % VELOCITY_U_3 ), &
             T    =  FV ( :, F % TEMPERATURE ), &
             N_0  =  Rho_0 / m_b, &
             T_0  =  T_0 )

    end associate !-- FV
    end associate !-- L, etc.

  end subroutine SetFluid


  subroutine SetRadiation ( MT, R )

    class ( MarshakThermalizationForm ), intent ( inout ) :: &
      MT
    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      R

    associate &
      (      L       =>  MT % BoxLength, &
           Rho_0     =>  MT % MassDensity, &
         Kappa       =>  MT % SpecificOpacity, &
        Lambda       =>  MT % MeanFreePath, &
           Tau       =>  MT % OpticalDepth, &
             t_Diff  =>  MT % DiffusionTime, &
             c       =>  CONSTANT % SPEED_OF_LIGHT )
    associate &
      ( RV  =>  R % Storage_GS % Value )

    Lambda  =  1.0_KDR  /  ( Rho_0 * Kappa )
    Tau     =  L / Lambda
    t_Diff  =  L * Tau / c

    call SetRadiationKernel &
           ( J      =  RV ( :, R % ENERGY_DENSITY_C ), &
             H_1    =  RV ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
             H_2    =  RV ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
             H_3    =  RV ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
             T_0    =  MT % TemperatureRadiation, &
             a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN )

    end associate !-- RV
    end associate !-- L, etc.

  end subroutine SetRadiation


  subroutine SetFluidKernel ( N, V_1, V_2, V_3, T, N_0, T_0 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      V_1, V_2, V_3, &
      T
    real ( KDR ), intent ( in ) :: &
      N_0, &
      T_0

    N    =  N_0
    V_1  =  0.0_KDR
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR
    T    =  T_0

  end subroutine SetFluidKernel


  subroutine SetRadiationKernel ( J, H_1, H_2, H_3, T_0, a )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3
    real ( KDR ), intent ( in ) :: &
      T_0, &
      a

    J    =  a  *  T_0 ** 4
    H_1  =  0.0_KDR
    H_2  =  0.0_KDR
    H_3  =  0.0_KDR

  end subroutine SetRadiationKernel


end module MarshakThermalization_Form
