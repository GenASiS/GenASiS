module FishboneMoncrief_Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  use Fluid_P_P__Form
  use Fluid_ASC__Form
  use Tally_FM__Form

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: FishboneMoncriefForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
  end type FishboneMoncriefForm

    private :: &
      SetFluid, &
      ApplySources, &
      ApplySourcesKernel

    real ( KDR ), private :: &
      AngularMomentumParameter, &
      CentralMass, &
      RadiusInner, &
      AdiabaticIndex, &
      DensityMax, &
      AtmosphereParameter

contains


  subroutine Initialize ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iB  !-- iBoundary
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      RadiusMin, &
      RadiusOuter, &
      SpecificAngularMomentum, &
      FinishTime
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      MassDensityUnit, &
      EnergyDensityUnit, &
      MassUnit, &
      EnergyUnit, &
      MomentumUnit, &
      AngularMomentumUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit, &
      VelocityUnit
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    associate &
      ( Kappa  => AngularMomentumParameter, &  !-- Between 1 and 2
        M      => CentralMass, &
        R_In   => RadiusInner, &
        Gamma  => AdiabaticIndex, &
        N_Max  => DensityMax, &
        AP     => AtmosphereParameter, &
        R_Out  => RadiusOuter, &
        L      => SpecificAngularMomentum, &
        G      => CONSTANT % GRAVITATIONAL, &
        c      => CONSTANT % SPEED_OF_LIGHT )

    Kappa  =  1.85_KDR
    M      =  3.0_KDR * UNIT % SOLAR_MASS
    call PROGRAM_HEADER % GetParameter ( Kappa, 'AngularMomentumParameter' )
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )

    R_In   =  6.0_KDR * G * M  /  c ** 2
    Gamma  =  1.4_KDR
    N_Max  =  1.0e12_KDR * UNIT % MASS_DENSITY_CGS
    AP     =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter ( R_In,  'RadiusInner' )
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( N_Max, 'DensityMax' )
    call PROGRAM_HEADER % GetParameter ( AP,    'AtmosphereParameter' )

    R_Out  =  R_In * Kappa / ( 2.0_KDR - Kappa )
    L      =  sqrt ( Kappa * G * M * R_In )

    call Show ( 'FishboneMoncrief parameters' )
    call Show ( Kappa, 'Kappa' ) 
    call Show ( M,     UNIT % SOLAR_MASS, 'M' )
    call Show ( R_In,  UNIT % KILOMETER, 'R_In' )
    call Show ( R_Out, UNIT % KILOMETER, 'R_Out' )
    call Show ( L,     UNIT % KILOMETER * UNIT % SPEED_OF_LIGHT, 'L' )
    call Show ( Gamma, 'Gamma' )
    call Show ( N_Max, UNIT % MASS_DENSITY_CGS, 'N_Max' )

    RadiusMin  =  ( 2.0_KDR / 3.0_KDR ) * R_In
    call PROGRAM_HEADER % GetParameter ( RadiusMin, 'RadiusMin' )
    end associate !-- Kappa, etc.

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: FM % PositionSpace )
    select type ( PS => FM % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    nCells = [ 128, 128, 1 ]
    call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )

    CoordinateSystem = 'SPHERICAL'
    CoordinateUnit   = [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

    call PS % SetBoundaryConditionsFace &
           ( [ 'OUTFLOW', 'INFLOW ' ], iDimension = 1 )
    call PS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [   RadiusMin, 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ RadiusOuter,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    associate &
      ( dTheta => ( MaxCoordinate ( 2 )  -  MinCoordinate ( 2 ) ) &
                  / nCells ( 2 ) )
    Spacing = [ 'PROPORTIONAL', 'EQUAL       ', 'EQUAL       ' ]
    Ratio   = [ dTheta, 0.0_KDR, 0.0_KDR ]
    end associate !-- dTheta

    call PS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             nCellsOption = nCells )

    !-- Geometry of PositionSpace

    allocate ( FM % Geometry_ASC )
    associate ( GA => FM % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Fluid

    TimeUnit = UNIT % SECOND

    VelocityUnit ( 1 ) =  CoordinateUnit ( 1 ) / TimeUnit 
    VelocityUnit ( 2 ) =  CoordinateUnit ( 2 ) / TimeUnit
    VelocityUnit ( 3 ) =  CoordinateUnit ( 3 ) / TimeUnit
    MassDensityUnit    =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit  =  UNIT % MASS_DENSITY_CGS  &
                          *  UNIT % SPEED_OF_LIGHT ** 2

    MassUnit             =  UNIT % SOLAR_MASS
    EnergyUnit           =  MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumUnit         =  MassUnit  *  UNIT % SPEED_OF_LIGHT
    AngularMomentumUnit  =  CoordinateUnit ( 1 ) &
                            *  MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    allocate ( Fluid_ASC_Form :: FM % Current_ASC )
    select type ( FA => FM % Current_ASC )
    class is ( Fluid_ASC_Form )

    allocate ( Tally_FM_Form :: FA % TallyInterior )
    allocate ( Tally_FM_Form :: FA % TallyTotal )
    allocate ( Tally_FM_Form :: FA % TallyChange )
    allocate ( FA % TallyBoundaryLocal ( PS % nBoundaries ) )
    allocate ( FA % TallyBoundaryGlobal ( PS % nBoundaries ) )
    do iB = 1, PS % nBoundaries 
      allocate ( Tally_FM_Form :: FA % TallyBoundaryLocal ( iB ) % Element )
      allocate ( Tally_FM_Form :: FA % TallyBoundaryGlobal ( iB ) % Element )
    end do !-- iB
       
    call FA % Initialize &
           ( PS, 'POLYTROPIC', VelocityUnitOption = VelocityUnit, &
             MassDensityUnitOption = MassDensityUnit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             MassUnitOption = MassUnit, EnergyUnitOption = EnergyUnit, &
             MomentumUnitOption = MomentumUnit, &
             AngularMomentumUnitOption = AngularMomentumUnit )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FM % Step )
    select type ( S => FM % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources
    end select !-- S

    !-- Set fluid and initialize Integrator template

    FinishTime = 1.0e-3_KDR * UNIT % SECOND

    call SetFluid ( FM )
    call FM % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, FinishTimeOption = FinishTime )

    !-- Cleanup

    end select !-- FA
    end select !-- PS

  end subroutine Initialize


  subroutine Finalize ( FM )

    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

    call FM % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetFluid ( FM )

    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

    real ( KDR ) :: &
      EnthalpyMax, &
      PolytropicParameter
    real ( KDR ), dimension ( : ), allocatable :: &
      Enthalpy
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_P_Form ), pointer :: &
      F

    integer ( KDI ) :: &
      iB  !-- iBoundary

    associate &
      ( Kappa => AngularMomentumParameter, &
        M     => CentralMass, &
        R_In  => RadiusInner, &
        Gamma => AdiabaticIndex, &
        N_Max => DensityMax, &
        AP    => AtmosphereParameter, &
        W_Max => EnthalpyMax, &
        K     => PolytropicParameter, &
        GC    => CONSTANT % GRAVITATIONAL, &
        c     => CONSTANT % SPEED_OF_LIGHT )

    select type ( FA => FM % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_P ( )

    !-- Fluid

    call F % SetAdiabaticIndex ( Gamma )

    W_Max = GC * M / ( c ** 2  *  R_In )  &
            *  ( 0.5_KDR * ( Kappa  +  1.0_KDR / Kappa )  -  1.0_KDR )

    K  =  ( Gamma - 1.0_KDR ) / Gamma  &
          *  W_Max  /  N_Max ** ( Gamma - 1.0_KDR )

    allocate ( Enthalpy ( F % nValues ) )

    select type ( PS => FM % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (     W => Enthalpy, &
            N => F % Value ( :, F % COMOVING_DENSITY ), &
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
            E => F % Value ( :, F % INTERNAL_ENERGY ), &
            R => G % Value ( :, G % CENTER ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER ( 2 ) ) )

    W = max ( 0.0_KDR, &
              GC * M / ( c ** 2  *  R_In )  &
              * ( R_In / R  -  1.0_KDR  +  0.5_KDR * Kappa &
                  -  0.5_KDR * Kappa  *  R_In ** 2  &
                     /  ( R * sin ( Theta ) ) ** 2 ) )

    N  =  ( ( Gamma - 1.0_KDR ) / ( Gamma * K ) * W ) &
          ** ( 1.0_KDR / ( Gamma - 1.0_KDR ) )

    E  =  K  *  N ** Gamma  /  ( Gamma - 1.0_KDR )

    V_1 = 0.0_KDR
    V_2 = 0.0_KDR

    where ( N > 0.0_KDR )
      V_3  =  sqrt ( Kappa * GC * M * R_In )  /  ( R * sin ( Theta ) ) ** 2
    elsewhere
      N    =    AP * N_Max  *  ( R / R_In ) ** ( - 1.5_KDR )
      V_1  =  - sqrt ( 2.0_KDR * GC * M / R )
      V_3  =    0.0_KDR
    end where

    call F % ComputeFromPrimitive ( G )

    !-- Tally

    select type ( TI => FA % TallyInterior )
    class is ( Tally_FM_Form )
       call TI % SetCentralMass ( M )
    end select !-- TI
    
    select type ( TT => FA % TallyTotal )
    class is ( Tally_FM_Form )
       call TT % SetCentralMass ( M )
    end select !-- TT
    
    select type ( TC => FA % TallyChange )
    class is ( Tally_FM_Form )
       call TC % SetCentralMass ( M )
    end select !-- TC
    
    do iB = 1, size ( FA % TallyBoundaryLocal )
      select type ( TB => FA % TallyBoundaryLocal ( iB ) % Element )
      class is ( Tally_FM_Form )
        call TB % SetCentralMass ( M )
      end select !-- TB
      select type ( TB => FA % TallyBoundaryGlobal ( iB ) % Element )
      class is ( Tally_FM_Form )
        call TB % SetCentralMass ( M )
      end select !-- TB
    end do !-- iB

    end associate !-- W, etc.
    end select !-- PS
    end select !-- FA
    end associate !-- Kappa, etc.
    nullify ( F, G )

  end subroutine SetFluid


  subroutine ApplySources ( S, Increment, Fluid, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iMomentum_1, &
      iEnergy
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call ApplySourcesCurvilinear_Fluid_P ( S, Increment, Fluid, TimeStep )

    select type ( F => Fluid )
    class is ( Fluid_P_Template )

    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    G => Chart % Geometry ( )

    call ApplySourcesKernel &
           ( Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iEnergy ), &
             Chart % IsProperCell, &
             F % Value ( :, F % COMOVING_DENSITY ), &
             F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             G % Value ( :, G % CENTER ( 1 ) ), &
             CONSTANT % GRAVITATIONAL, CentralMass, TimeStep ) 

    end select !-- Grid
    end select !-- F

    nullify ( G )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplySources

  
  subroutine ApplySourcesKernel &
               ( KV_M_1, KV_E, IsProperCell, N, V_1, R, G, M, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_M_1, &
      KV_E   
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, &
      R
    real ( KDR ) :: &
      G, &
      M, &
      dT

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ), dimension ( size ( N ) ) :: &
      F_1

    nV = size ( N )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      F_1 ( iV )  =  - G * M * N ( iV )  /  R ( iV ) ** 2
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      KV_M_1 ( iV )  =  KV_M_1 ( iV )  +  dT * F_1 ( iV )
      KV_E   ( iV )  =  KV_E   ( iV )  +  dT * F_1 ( iV ) * V_1 ( iV ) 
    end do
    !$OMP end parallel do

  end subroutine ApplySourcesKernel


end module FishboneMoncrief_Form
