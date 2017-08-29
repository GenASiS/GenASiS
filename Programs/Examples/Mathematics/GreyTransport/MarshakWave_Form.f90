module MarshakWave_Form

  !-- Vaytet et al. 2011

  use Basics
  use Mathematics
  use Fluid_P_NR__Form
  use Fluid_ASC__Form
  use RadiationMoments_Form
  use PhotonMoments_G__Form
  use RadiationMoments_ASC__Form
  use Interactions_Template
  use Interactions_MWV_1_G__Form
  use Interactions_MWV_2_G__Form
  use Interactions_MWV_3_G__Form
  use Interactions_ASC__Form
  use ApplyRelaxation_RM__Command

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_Template ) :: MarshakWaveForm
    integer ( KDI ) :: &
      RADIATION = 1, &
      FLUID     = 2         
    type ( Interactions_ASC_Form ), allocatable :: &
      Interactions_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type MarshakWaveForm

    private :: &
      SetRadiation, &
      SetFluid, &
      SetInteractions, &
      ApplySources_Radiation, &
      ApplySources_Fluid, &
      ComputeFluidSource_Radiation

      private :: &
        ApplySources_Fluid_Kernel, &
        ComputeFluidSource_Radiation_Kernel

    real ( KDR ), private :: &
      AdiabaticIndex, &
      SpecificHeatCapacity, &
      MeanMolecularWeight, &
      MassDensity, &
      Temperature, &
      TemperatureInner, &
      SpecificOpacity, &
      EnergyMax, &
      SoundSpeed
    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDF ) :: &
      InteractionsType
    type ( VariableGroupForm ), private, allocatable :: &
      FluidSource_Radiation
    class ( Fluid_P_NR_Form ), pointer :: &
      Fluid => null ( )
    class ( InteractionsTemplate ), pointer :: &
      Interactions => null ( )

contains


  subroutine Initialize ( MW, Name )

    class ( MarshakWaveForm ), intent ( inout ) :: &
      MW
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    real ( KDR ) :: &
      BoxLength, &
      MeanFreePath, &
      OpticalDepth, &
      DynamicalTime, &
      DiffusionTime, &
      FinishTime
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

    call PROGRAM_HEADER % GetParameter ( L,     'BoxLength' )
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( C_V,   'SpecificHeatCapacity' )
    call PROGRAM_HEADER % GetParameter ( N_0,   'MassDensity' )
    call PROGRAM_HEADER % GetParameter ( T_0,   'Temperature' )
    call PROGRAM_HEADER % GetParameter ( T_I,   'TemperatureInner' )
    call PROGRAM_HEADER % GetParameter ( Kappa, 'SpecificOpacity' )

    InteractionsType = 'MARSHAK_WAVE_VAYTET_1_GREY'
    call PROGRAM_HEADER % GetParameter ( InteractionsType, 'InteractionsType' )

    select case ( trim ( InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_2_GREY', 'MARSHAK_WAVE_VAYTET_3_GREY' )
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

    !-- Prepare for Currents

    MW % N_CURRENTS_PS = 2
    allocate ( MW % Current_ASC_1D ( MW % N_CURRENTS_PS ) )
    allocate ( MW % TimeStepLabel ( MW % N_CURRENTS_PS ) )
    MW % TimeStepLabel ( MW % RADIATION ) = 'Radiation'
    MW % TimeStepLabel ( MW % FLUID )     = 'Fluid'

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
    
    !-- RadiationMoments

    allocate &
      ( RadiationMoments_ASC_Form :: &
          MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    select type ( RA => MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA % Initialize &
           ( PS, 'PHOTONS_GREY', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             TemperatureUnitOption = TemperatureUnit )

    !-- Fluid

    allocate &
      ( Fluid_ASC_Form :: MW % Current_ASC_1D ( MW % FLUID ) % Element )
    select type ( FA => MW % Current_ASC_1D ( MW % FLUID ) % Element )
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

    allocate ( MW % Interactions_ASC )
    associate ( IA => MW % Interactions_ASC )
    call IA % Initialize &
           ( PS, InteractionsType, &
             LengthUnitOption = CoordinateUnit ( 1 ), &
             EnergyDensityUnitOption = EnergyDensityUnit )
    call RA % SetInteractions ( IA )

    !-- Step

    allocate ( Step_RK2_C_ASC_1D_Form :: MW % Step )
    select type ( S => MW % Step )
    class is ( Step_RK2_C_ASC_1D_Form )

    call S % Initialize ( MW % Current_ASC_1D, Name )

    S % ApplySources_1D ( MW % RADIATION ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplySources_1D ( MW % FLUID ) % Pointer &
      =>  ApplySources_Fluid
    S % ApplyRelaxation_1D ( MW % RADIATION ) % Pointer &
      =>  ApplyRelaxation_RM
    S % HarvestIncrement_1D ( MW % RADIATION ) % Pointer &
      =>  ComputeFluidSource_Radiation

    end select !-- S

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
        c      => CONSTANT % SPEED_OF_LIGHT )

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

    select case ( trim ( InteractionsType ) )
    case ( 'MARSHAK_WAVE_VAYTET_2_GREY', 'MARSHAK_WAVE_VAYTET_3_GREY' )
      call Show ( EnergyMax, UNIT % ELECTRON_VOLT, 'EnergyMax' )
    end select !-- InteractionsType

    end associate !-- c_s, etc.

    !-- Initialize template

    FinishTime  =  1.36e-7_KDR  *  UNIT % SECOND

    call MW % InitializeTemplate_C_1D_PS &
           ( Name, TimeUnitOption = TimeUnit, FinishTimeOption = FinishTime )

    !-- Cleanup

    end associate !-- IA
    end select !-- FA
    end select !-- RA
    end select !-- PS
    end associate !-- L, etc.

  end subroutine Initialize


  impure elemental subroutine Finalize ( MW )
    
    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    if ( allocated ( MW % Interactions_ASC ) ) &
      deallocate ( MW % Interactions_ASC )

    call MW % FinalizeTemplate_C_1D_PS ( )

  end subroutine Finalize


  subroutine SetRadiation ( MW )

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( PhotonMoments_G_Form ), pointer :: &
      R

    associate &
      ( T_0 => Temperature, &
        T_I => TemperatureInner )

    select type ( RA => MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    class is ( RadiationMoments_ASC_Form )
    R => RA % PhotonMoments_G ( )

    select type ( PS => MW % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (   J => R % Value ( :, R % COMOVING_ENERGY ), &
        H_1 => R % Value ( :, R % COMOVING_MOMENTUM_U ( 1 ) ), &
        H_2 => R % Value ( :, R % COMOVING_MOMENTUM_U ( 2 ) ), &
        H_3 => R % Value ( :, R % COMOVING_MOMENTUM_U ( 3 ) ), &
          X => G % Value ( :, G % CENTER ( 1 ) ), &
          Y => G % Value ( :, G % CENTER ( 2 ) ), &
          Z => G % Value ( :, G % CENTER ( 3 ) ), &
          a => 4.0_KDR * CONSTANT % STEFAN_BOLTZMANN )

    J  =  a  *  T_0 ** 4

    H_1  =  0.0_KDR
    H_2  =  0.0_KDR
    H_3  =  0.0_KDR

    where ( X < MinCoordinate ( 1 ) .or. Y < MinCoordinate ( 2 ) &
            .or. Z < MinCoordinate ( 3 ) )
      J  =  a  *  T_I ** 4
    end where

    call R % ComputeFromPrimitive ( G )

    end associate !-- J, etc.
    end select !-- PS
    end select !-- RA
    end associate !-- T_0
    nullify ( R, G )

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

    select type ( FA => MW % Current_ASC_1D ( MW % FLUID ) % Element )
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

    type ( MarshakWaveForm ), intent ( inout ) :: &
      MW

    class ( Fluid_P_NR_Form ), pointer :: &
      F
    class ( PhotonMoments_G_Form ), pointer :: &
      R
    class ( InteractionsTemplate), pointer :: &
      I

    associate ( IA => MW % Interactions_ASC )

    select type ( RA => MW % Current_ASC_1D ( MW % RADIATION ) % Element )
    class is ( RadiationMoments_ASC_Form )
    R => RA % PhotonMoments_G ( )

    select type ( FA => MW % Current_ASC_1D ( MW % FLUID ) % Element )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_NR ( )

    I => IA % Interactions ( )
    select type ( I )
    type is ( Interactions_MWV_1_G_Form )
      call I % Set ( R, F, SpecificOpacity )
    type is ( Interactions_MWV_2_G_Form )
      call I % Set ( R, F, SpecificOpacity, EnergyMax )
    type is ( Interactions_MWV_3_G_Form )
      call I % Set ( R, F, SpecificOpacity, EnergyMax, Temperature )
    class default
      call Show ( 'Interactions type not recognized', CONSOLE % ERROR )
      call Show ( 'MarshakWave_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetInteractions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I

    call I % Compute ( R )

    !-- Module variable for accessibility in ApplySources_Fluid below
    Interactions => I

    !-- Storage for fluid source terms

    allocate ( FluidSource_Radiation )
    call FluidSource_Radiation % Initialize &
           ( [ F % nValues, F % N_CONSERVED ], ClearOption = .true. )

    end select !-- FA
    end select !-- RA
    end associate !-- IA
    nullify ( F, R, I )

  end subroutine SetInteractions


  subroutine ApplySources_Radiation &
               ( S, Sources_RM, Increment, Radiation, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    !-- No sources applied here; just an occasion to compute interactions
    !   to be used in relaxation.

    call Interactions % Compute ( Radiation )

  end subroutine ApplySources_Radiation

  
  subroutine ApplySources_Fluid &
               ( S, Sources_RM, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
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
             FluidSource_Radiation % Value ( :, iEnergy ), &
             FluidSource_Radiation % Value ( :, iMomentum_1 ) )

    call Clear ( FluidSource_Radiation % Value )

    end select !-- Chart
    end select !-- F

  end subroutine ApplySources_Fluid

  
  subroutine ComputeFluidSource_Radiation ( S, Increment, Radiation, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( inout ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iEnergy_F, &
      iEnergy_R, &
      iMomentum_1_F, &
      iMomentum_1_R

    select type ( R => Radiation )
    class is ( PhotonMoments_G_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    associate &
      ( I => Interactions, &
        F => Fluid )

    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy_F )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), &
                  iMomentum_1_F )
    call Search ( R % iaConserved, R % CONSERVED_ENERGY, iEnergy_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_D ( 1 ), &
                  iMomentum_1_R )

    !-- Taking shortcuts on conserved vs. comoving here
    call ComputeFluidSource_Radiation_Kernel &
           ( FluidSource_Radiation % Value ( :, iEnergy_F ), &
             FluidSource_Radiation % Value ( :, iMomentum_1_F ), &
             Chart % IsProperCell, &
             I % Value ( :, I % EMISSIVITY_H ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ), &
             R % Value ( :, R % COMOVING_ENERGY ), &
             Increment % Value ( :, iEnergy_R ), &
             R % Value ( :, R % CONSERVED_MOMENTUM_D ( 1 ) ), &
             Increment % Value ( :, iMomentum_1_R ), &
             CONSTANT % SPEED_OF_LIGHT, TimeStep ) 

    end associate !-- I, etc.
    end select !-- Chart
    end select !-- R

  end subroutine ComputeFluidSource_Radiation


  subroutine ApplySources_Fluid_Kernel &
               ( K_G, K_S_1, IsProperCell, FS_R_G, FS_R_S_1 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_G, &
      K_S_1
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FS_R_G, &
      FS_R_S_1

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( K_G )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      K_G ( iV )    =  K_G ( iV )    +  FS_R_G ( iV )
      K_S_1 ( iV )  =  K_S_1 ( iV )  +  FS_R_S_1 ( iV )
    end do
    !$OMP end parallel do

  end subroutine ApplySources_Fluid_Kernel


  subroutine ComputeFluidSource_Radiation_Kernel &
               ( FS_R_G, FS_R_S_1, IsProperCell, Xi_J, Chi_J, Chi_H, &
                 J, dJ, H, dH, c, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_G, &
      FS_R_S_1
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H, &
      J,  &
      dJ, &
      H,  &
      dH
    real ( KDR ) :: &
      c, &
      dT

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( FS_R_G )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      FS_R_G ( iV )  &
        =  FS_R_G ( iV )  &
           -  c * dT  *  ( Xi_J ( iV )  &
                           -  Chi_J ( iV ) * ( J ( iV ) + dJ ( iV ) ) ) 
      FS_R_S_1 ( iV )  &
        =  FS_R_S_1 ( iV )  &
           +  c * dT  *  Chi_H ( iV )  *  ( H ( iV ) + dH ( iV ) )
    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_Radiation_Kernel


end module MarshakWave_Form
