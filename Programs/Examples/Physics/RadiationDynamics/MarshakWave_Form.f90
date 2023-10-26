module MarshakWave_Form

  !-- Vaytet et al. 2011

  use GenASiS
  use Interactions_MWV_1__Form

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: MarshakWaveForm
    real ( KDR ) :: &
      BoxLength, &
      AdiabaticIndex, &
      SpecificHeatCapacity, &  !-- per unit mass
      MassDensity, &
      Temperature, &
      TemperatureInner, &
      SpecificOpacity, &
      SpecificOpacityFloor!, &
    !   EnergyMax, &
    !   SoundSpeed
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
    procedure, public, pass :: &
      ShowParameters
  end type MarshakWaveForm

    private :: &
      InitializeUniverse, &
      SetInitial

    !   private :: &
    !     SetFluid, &
    !     SetRadiation

    !     private :: &
    !       SetFluidKernel, &
    !       SetPerturbationKernel


contains


  subroutine Initialize_MW ( MW, FormalismType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( MW % Type  ==  '' ) &
      MW % Type  =  'a Thermalization'

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
    call Show ( U % SpecificOpacityFloor, &
                UNIT % CENTIMETER ** 2 / UNIT % GRAM, &
                'SpecificOpacityFloor' )
                
  end subroutine ShowParameters


  subroutine InitializeUniverse ( MW, FormalismType, Name )

    class ( MarshakWaveForm ), intent ( inout ), target :: &
      MW
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    integer ( KDI ) :: &
      iD

    !-- Units

    allocate ( MW % Units_F ( 1 ) )
    call MW % Units_F ( 1 ) % Initialize ( TypeOption = 'CGS' )

    !-- Position space parameters

    MW % BoxLength  =  25.0_KDR  *  UNIT % CENTIMETER
    call PROGRAM_HEADER % GetParameter ( MW % BoxLength, 'BoxLength' )

    MW % MinCoordinate  =  0.0_KDR
    MW % MaxCoordinate  =  MW % BoxLength

    !-- Interactions

    allocate ( Interactions_MWV_1_Form :: MW % Interactions_BM )

    !-- Initialization

    call MW % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'PHOTONS' ], &
             FormalismType = FormalismType, &
             Name = Name, &
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
             ( [ 'INFLOW', 'INFLOW' ], iC = 1, iD = iD )
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
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )

    !-- FinishTime

    I % T_Finish  =  1.36e-7_KDR  *  UNIT % SECOND

    !-- Parameters

    associate &
      ( Gamma      =>  MW % AdiabaticIndex, &
        C_V        =>  MW % SpecificHeatCapacity, &
        Rho_0      =>  MW % MassDensity, &
        T_0        =>  MW % Temperature, &
        T_I        =>  MW % TemperatureInner, &
        Kappa      =>  MW % SpecificOpacity, &
        Kappa_Min  =>  MW % SpecificOpacityFloor )

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

    end associate !-- Gamma, etc.

    end select !-- F
    end select !-- I
    end select !-- T

  end subroutine SetInitial


end module MarshakWave_Form
