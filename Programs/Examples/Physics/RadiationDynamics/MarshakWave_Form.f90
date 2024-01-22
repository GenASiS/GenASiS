module MarshakWave_Form

  !-- Vaytet et al. 2011

  use GenASiS
  use Interactions_MWV_1__Form
  use Interactions_MWV_2__Form
  use Interactions_MWV_3__Form

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: MarshakWaveForm
    real ( KDR ) :: &  !-- Primary parameters
      BoxLength, &
      AdiabaticIndex, &
      SpecificHeatCapacity, &  !-- per unit mass
      MassDensity, &
      Temperature, &
      TemperatureInner, &
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
      Initialize_MW
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type MarshakWaveForm

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


  subroutine Initialize_MW ( MW, FormalismType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( MW % Type  ==  '' ) &
      MW % Type  =  'a MarshakWave'

    call InitializeUniverse ( MW, FormalismType, Name )
 
  end subroutine Initialize_MW


  impure elemental subroutine Finalize ( MW )

    type ( MarshakWaveForm ), intent ( inout ), target :: &
      MW

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( MarshakWaveForm ), intent ( in ) :: &
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
    call Show ( U % Temperature, &
                UNIT % KELVIN, &
                'Temperature' )
    call Show ( U % TemperatureInner, &
                UNIT % KELVIN, &
                'TemperatureInner' )
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


  subroutine InitializeUniverse ( MW, FormalismType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    integer ( KDI ) :: &
      iD

    !-- Position space parameters

    MW % BoxLength  =  25.0_KDR  *  UNIT % CENTIMETER
    call PROGRAM_HEADER % GetParameter ( MW % BoxLength, 'BoxLength' )

    MW % MinCoordinate  =  0.0_KDR
    MW % MaxCoordinate  =  MW % BoxLength

    ! !-- Momentum space parameters

    ! !-- Geometric spacing
    ! MinWidthEnergy  =  3.0e1_KDR  *  UNIT % KELVIN
    ! MaxEnergy       =  1.0e4_KDR  *  UNIT % KELVIN

    ! !-- Compactified spacing
    ! EnergyScale     =  1.0e3_KDR  *  UNIT % KELVIN


    !-- Interactions

    MW % InteractionsType = 'MARSHAK_WAVE_VAYTET_1'
    call PROGRAM_HEADER % GetParameter &
           ( MW % InteractionsType, 'InteractionsType' )

    select case ( trim ( MW % InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_1' )
      allocate ( Interactions_MWV_1_Form :: MW % Interactions_BM )
    case ( 'MARSHAK_WAVE_VAYTET_2' )
      allocate ( Interactions_MWV_2_Form :: MW % Interactions_BM )
    case ( 'MARSHAK_WAVE_VAYTET_3' )
      allocate ( Interactions_MWV_3_Form :: MW % Interactions_BM )
    case default
      call Show ( 'InteractionsType not recognized', CONSOLE % ERROR )
      call Show ( 'MarshakWave_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeUniverse', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    !-- Initialization

    call MW % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'PHOTONS' ], &
             FormalismType = FormalismType, &
             Name = Name, &
             UnitsTypeOption = 'CGS', &
             MinCoordinateOption = MW % MinCoordinate, &
             MaxCoordinateOption = MW % MaxCoordinate, &
             nCellsPositionOption = [ 128, 128, 128 ] )

    !-- Boundary conditions

    select type ( I  =>  MW % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    select type ( I  =>  MW % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    associate &
      ( R  =>  I % CurrentSet_X_1D )
    do iD  =  1, 3
      call R % SetBoundaryConditionsFace &
             ( [ 'INFLOW', 'INFLOW' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- R
    end select !-- I

    !-- Pointers

    MW % Integrator % SetInitial    =>  SetInitial
    MW % Integrator % System        =>  MW
    
  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( MW  =>  I % System )
      class is ( MarshakWaveForm )
    select type ( I )
      class is ( Integrator_CS_1D_CS_Form )

    !-- FinishTime

!    I % T_Finish  =  1.36e-7_KDR  *  UNIT % SECOND
!-- More diffusive
    I % T_Finish  =  1.36e-6_KDR  *  UNIT % SECOND

    !-- Parameters

    associate &
      ( Gamma      =>  MW % AdiabaticIndex, &
        C_V        =>  MW % SpecificHeatCapacity, &
        Rho_0      =>  MW % MassDensity, &
        T_0        =>  MW % Temperature, &
        T_I        =>  MW % TemperatureInner, &
        Kappa      =>  MW % SpecificOpacity, &
        Kappa_Min  =>  MW % SpecificOpacityMin, &
        E_Max      =>  MW % EnergyMax )

    Gamma      =  1.4_KDR
    C_V        =  1.0_KDR     *  UNIT % ERG / UNIT % KELVIN / UNIT % GRAM
    Rho_0      =  1.0e-3_KDR  *  UNIT % MASS_DENSITY_CGS
    T_0        =  3.0e2_KDR   *  UNIT % KELVIN
    T_I        =  1.0e3_KDR   *  UNIT % KELVIN
!    Kappa      =  1.0e3_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
!-- More diffusive
    Kappa      =  1.0e4_KDR   *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    Kappa_Min  =  10.0_KDR    *  UNIT % CENTIMETER ** 2 / UNIT % GRAM
    E_Max      =  0.620_KDR   *  UNIT % ELECTRON_VOLT

    call PROGRAM_HEADER % GetParameter ( Gamma,     'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,       'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( Rho_0,     'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,       'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,       'TemperatureInner' )
    call PROGRAM_HEADER % GetParameter ( Kappa,     'SpecificOpacity' )
    call PROGRAM_HEADER % GetParameter ( Kappa_Min, 'SpecificOpacityMin' )
    call PROGRAM_HEADER % GetParameter ( E_Max,     'EnergyMax' )

    end associate !-- Gamma, etc.

    !-- Fluid

    select type ( F  =>  I % CurrentSet_X )
    class is ( Fluid_P_I_Form )
      call SetFluid ( MW, F )
      call F % SetUseInitialTemperature ( .true. )
    end select !-- F

    !-- Radiation

    select type ( I )
    class is ( Integrator_CS_1D_BM_CS_Form )

      select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )

      call SetRadiation ( MW, R )

      end select !-- R

    class default
      call Show ( 'Integrator type not recognized', CONSOLE % ERROR )
      call Show ( 'MarshakWaveForm', 'module', CONSOLE % ERROR )
      call Show ( 'SetInitial', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I    

    !-- Interactions

    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( PhotonMoments_G_Form )

    select type ( Intrctns  =>  MW % Interactions_BM )
    class is ( Interactions_MWV_3_Form )
       call Intrctns % SetSpecificOpacity ( MW % SpecificOpacity )
       call Intrctns % SetEnergyMax ( MW % EnergyMax )
       call Intrctns % SetTemperatureScale ( MW % Temperature )
    class is ( Interactions_MWV_2_Form )
       call Intrctns % SetSpecificOpacity ( MW % SpecificOpacity )
       call Intrctns % SetEnergyMax ( MW % EnergyMax )
    class is ( Interactions_MWV_1_Form )
       call Intrctns % SetSpecificOpacity ( MW % SpecificOpacity )
    end select !-- Intrctns
    
    end select !-- R
    end select !-- I

    !-- Cleanup

    end select !-- I
    end select !-- MW

  end subroutine SetInitial


  subroutine SetFluid ( MW, F )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    real ( KDR ) :: &
      N_0, &
      E_0, &
      P_0

    associate &
      (     L      =>  MW % BoxLength, &
        Gamma      =>  MW % AdiabaticIndex, &
            C_V    =>  MW % SpecificHeatCapacity, &
          Rho_0    =>  MW % MassDensity, &
            T_0    =>  MW % Temperature, &
            c_s    =>  MW % SoundSpeed, &
            t_Dyn  =>  MW % DynamicalTime, &
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


  subroutine SetRadiation ( MW, R )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      R

    associate &
      ( G  =>  R % Geometry )
    associate &
      ( RV  =>  R % Storage_GS % Value, & 
        GV  =>  G % Storage_GS % Value )
    associate &
      (      L       =>  MW % BoxLength, &
           Rho_0     =>  MW % MassDensity, &
         Kappa       =>  MW % SpecificOpacity, &
        Lambda       =>  MW % MeanFreePath, &
           Tau       =>  MW % OpticalDepth, &
             t_Diff  =>  MW % DiffusionTime, &
             c       =>  CONSTANT % SPEED_OF_LIGHT )

    Lambda  =  1.0_KDR  /  ( Rho_0 * Kappa )
    Tau     =  L / Lambda
    t_Diff  =  L * Tau / c

    call SetRadiationKernel &
           ( J      =  RV ( :, R % ENERGY_DENSITY_C ), &
             H_1    =  RV ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
             H_2    =  RV ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
             H_3    =  RV ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
             X_1    =  GV ( :, G % CENTER_U_1 ), & 
             X_2    =  GV ( :, G % CENTER_U_2 ), & 
             X_3    =  GV ( :, G % CENTER_U_3 ), &
             X_Min  =  MW % MinCoordinate, &
             T_0    =  MW % Temperature, &
             T_I    =  MW % TemperatureInner, &
             a      =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN )

    end associate !-- L, etc.
    end associate !-- RV, etc.
    end associate !-- G

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


  subroutine SetRadiationKernel &
               ( J, H_1, H_2, H_3, X_1, X_2, X_3, X_Min, T_0, T_I, a )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1, X_2, X_3
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      X_Min
    real ( KDR ), intent ( in ) :: &
      T_0, T_I, &
      a

      where (      X_1  <  X_Min ( 1 )  &
             .or.  X_2  <  X_Min ( 2 )  &
             .or.  X_3  <  X_Min ( 3 ) )

        J  =  a  *  T_I ** 4

      elsewhere

        J  =  a  *  T_0 ** 4

      end where
    
  end subroutine SetRadiationKernel


end module MarshakWave_Form
