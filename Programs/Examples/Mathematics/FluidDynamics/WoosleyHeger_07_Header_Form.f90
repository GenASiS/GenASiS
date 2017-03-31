module WoosleyHeger_07_Header_Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  use Fluid_P_MHN__Form
  use Fluid_ASC__Form
  use Tally_F_P_MHN__Form

  implicit none
  private

  type, public :: WoosleyHeger_07_HeaderForm
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
      VelocityUnit
    class ( IntegratorTemplate ), pointer :: &
      Integrator => null ( )
    class ( Fluid_ASC_Form ), pointer :: &
      Fluid_ASC => null ( )
  contains
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SetFluid
  end type WoosleyHeger_07_HeaderForm

    private :: &
      PrepareInterpolation

    integer ( KDI ), private, parameter :: &
      iRADIUS_TS            = 2, &  !-- must match the profile file columns
      iRADIAL_VELOCITY_TS   = 3, &
      iDENSITY_TS           = 4, &
      iTEMPERATURE_TS       = 5, &
      iSPECIFIC_ENERGY_TS   = 10, &
      iELECTRON_FRACTION_TS = 11
    integer ( KDI ), private, parameter :: &
      iRADIAL_VELOCITY_SI   = 1, &  !-- spline interpolation
      iDENSITY_SI           = 2, &
      iTEMPERATURE_SI       = 3, &
      iSPECIFIC_ENERGY_SI   = 4, &
      iELECTRON_FRACTION_SI = 5

  public :: &
    ApplySourcesGravity

    private :: &
      ComputeGravity, &
      ApplySourcesKernel

      private :: &
        ComputeGravityEdge, &
        ComputeGravityCenter!, &
  !     LocalMax

    real ( KDR ), dimension ( : ), allocatable, private :: &
      Potential, &    !-- Edge
      Potential_C, &  !-- Center
      Acceleration_C  !-- Center

contains


  subroutine InitializePositionSpace ( WHH, I )

    class ( WoosleyHeger_07_HeaderForm ), intent ( inout ) :: &
      WHH
    class ( IntegratorTemplate ), intent ( in ), target :: &
      I

    integer ( KDI ) :: &
      nCellsCore, &
      nCellsRadius, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    WHH % Integrator => I

    allocate ( Atlas_SC_Form :: WHH % Integrator % PositionSpace )
    select type ( PS => WHH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    CoordinateSystem = 'SPHERICAL'
    WHH % CoordinateUnit   = [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

    call PS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
    if ( PS % nDimensions > 1 ) &
      call PS % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    if ( PS % nDimensions > 2 ) &
      call PS % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 3 )

    RadiusMax = 1.0e4_KDR * UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'RadiusMax' )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [   0.0_KDR, 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    RadiusCore = 40.0_KDR  *  UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( RadiusCore, 'RadiusCore' )

    nCellsCore = 32  !-- Number of central cells with equal spacing
    call PROGRAM_HEADER % GetParameter ( nCellsCore, 'nCellsCore' )

    call Show ( 'Mesh core parameters' )
    call Show ( RadiusCore, UNIT % KILOMETER, 'RadiusCore' )
    call Show ( nCellsCore, 'nCellsCore' )
    call Show ( RadiusCore / nCellsCore, UNIT % KILOMETER, 'CellWidthCore' )

    nCellsRadius = 6.5 * nCellsCore  !-- Aiming for roughly 10,000 km
    call PROGRAM_HEADER % GetParameter ( nCellsRadius, 'nCellsRadius' )

    nCellsPolar     = 3 * nCellsCore
    nCellsAzimuthal = 2 * nCellsPolar

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  RadiusCore

    nCells = [ nCellsRadius, 1, 1 ]
    if ( PS % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( PS % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    call PS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = WHH % CoordinateUnit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nEqualOption = nCellsCore )

    !-- Geometry of PositionSpace

    allocate ( WHH % Integrator % Geometry_ASC )
    associate ( GA => WHH % Integrator % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    WHH % TimeUnit = UNIT % SECOND

    end select !-- PS

  end subroutine InitializePositionSpace


  subroutine InitializeFluid ( WHH, FA )

    class ( WoosleyHeger_07_HeaderForm ), intent ( inout ) :: &
      WHH
    type ( Fluid_ASC_Form ), intent ( inout ), target :: &
      FA

    WHH % Fluid_ASC => FA

    WHH % VelocityUnit       =  WHH % CoordinateUnit  /  WHH % TimeUnit 
    WHH % MassDensityUnit    =  UNIT % MASS_DENSITY_CGS
    WHH % EnergyDensityUnit  =  UNIT % MASS_DENSITY_CGS  &
                                *  UNIT % SPEED_OF_LIGHT ** 2
    WHH % TemperatureUnit    =  UNIT % MEV

    WHH % MassUnit             =  UNIT % SOLAR_MASS
    WHH % EnergyUnit           =  WHH % MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    WHH % MomentumUnit         =  WHH % MassUnit  *  UNIT % SPEED_OF_LIGHT
    WHH % AngularMomentumUnit  =  WHH % CoordinateUnit ( 1 ) &
                                  *  WHH % MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    select type ( PS => WHH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    call FA % Initialize &
           ( PS, 'MEAN_HEAVY_NUCLEUS', &
             VelocityUnitOption        = WHH % VelocityUnit, &
             MassDensityUnitOption     = WHH % MassDensityUnit, &
             EnergyDensityUnitOption   = WHH % EnergyDensityUnit, &
             TemperatureUnitOption     = WHH % TemperatureUnit, &
             MassUnitOption            = WHH % MassUnit, &
             EnergyUnitOption          = WHH % EnergyUnit, &
             MomentumUnitOption        = WHH % MomentumUnit, &
             AngularMomentumUnitOption = WHH % AngularMomentumUnit, &
             ShockThresholdOption = 1.0_KDR )

    end select !-- PS

  end subroutine InitializeFluid


  subroutine SetFluid ( WHH )

    class ( WoosleyHeger_07_HeaderForm ), intent ( inout ) :: &
      WHH

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iB     !-- iBoundary
    real ( KDR ) :: &
      SE  !-- SpecificEnergy
    type ( SplineInterpolationForm ), dimension ( 5 ) :: &
      SI
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_MHN_Form ), pointer :: &
      F

    call PrepareInterpolation ( SI )

    associate ( FA => WHH % Fluid_ASC )
    F => FA % Fluid_P_MHN ( )

    select type ( PS => WHH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (     N => F % Value ( :, F % COMOVING_DENSITY ), &
            T => F % Value ( :, F % TEMPERATURE ), &
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
            E => F % Value ( :, F % INTERNAL_ENERGY ), &
          Y_E => F % Value ( :, F % ELECTRON_FRACTION ), &
            R => G % Value ( :, G % CENTER ( 1 ) ) )

    do iV = 1, size ( N )

      call SI ( iDENSITY_SI ) % Evaluate ( R ( iV ), N ( iV ) ) 
      call SI ( iTEMPERATURE_SI ) % Evaluate ( R ( iV ) , T ( iV ) ) 
      call SI ( iRADIAL_VELOCITY_SI ) % Evaluate ( R ( iV ), V_1 ( iV ) ) 
      call SI ( iSPECIFIC_ENERGY_SI ) % Evaluate ( R ( iV ), SE ) 
      call SI ( iELECTRON_FRACTION_SI ) % Evaluate ( R ( iV ), Y_E ( iV ) ) 

      E ( iV )  =  SE  *  N ( iV )

    end do !-- iV

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    call F % ComputeFromPrimitive ( G )

    !-- Gravitational potential storage

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )
      allocate ( Potential ( 0 : C % nCells ( 1 ) ) )  !-- global edges
      allocate ( Potential_C ( G % nValues ) )  !-- local centers
      allocate ( Acceleration_C ( G % nValues ) )  !-- local centers
    end select !-- C

    !-- Tally

    select type ( TI => FA % TallyInterior )
    class is ( Tally_F_P_MHN_Form )
      call TI % SetGravitationalPotential ( Potential_C )
      TI % ComputeGravitationalPotential => ComputeGravity
    end select !-- TI
    
    select type ( TT => FA % TallyTotal )
    class is ( Tally_F_P_MHN_Form )
      call TT % SetGravitationalPotential ( Potential_C )
    end select !-- TT
    
    select type ( TC => FA % TallyChange )
    class is ( Tally_F_P_MHN_Form )
      call TC % SetGravitationalPotential ( Potential_C )
    end select !-- TC
    
    do iB = 1, size ( FA % TallyBoundaryLocal )
      select type ( TB => FA % TallyBoundaryLocal ( iB ) % Element )
      class is ( Tally_F_P_MHN_Form )
        call TB % SetGravitationalPotential ( Potential_C )
      end select !-- TB
      select type ( TB => FA % TallyBoundaryGlobal ( iB ) % Element )
      class is ( Tally_F_P_MHN_Form )
        call TB % SetGravitationalPotential ( Potential_C )
      end select !-- TB
    end do !-- iB

    end associate !-- N, etc.
    end select !-- PS
    end associate !-- FA
    nullify ( F, G )

  end subroutine SetFluid


  subroutine PrepareInterpolation ( SI )

    type ( SplineInterpolationForm ), dimension ( 5 ), intent ( inout ) :: &
      SI

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Slope_N, &
      Slope_T, &
      Slope_SE, &
      Slope_Y_E
    real ( KDR ), dimension ( : ), allocatable :: &
      RC, &   !-- RadiusCenter
      dRC, &  !-- WidthCenter
      Radius, &
      RadialVelocity, &
      Density, &
      Temperature, &
      SpecificEnergy, &
      ElectronFraction
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Profile
    character ( LDF ) :: &
      Path, &
      Filename
    type ( TableStreamForm ) :: &
      TS

    Path = '../Parameters/'
    Filename = 'WH07_S12_08.d.stripped'

    call TS % Initialize &
           ( Filename, PROGRAM_HEADER % Communicator % Rank, &
             PathOption = Path )
    call TS % Read ( Profile, oRowOption = 2 )

    !-- Set "edge" values

    associate &
      (   R => Profile ( :, iRADIUS_TS ), &             !-- cell outer edge
        V_R => Profile ( :, iRADIAL_VELOCITY_TS ), &    !-- cell outer edge
          N => Profile ( :, iDENSITY_TS ), &            !-- cell center
          T => Profile ( :, iTEMPERATURE_TS ), &        !-- cell center
        Y_E => Profile ( :, iELECTRON_FRACTION_TS ), &  !-- cell center
         SE => Profile ( :, iSPECIFIC_ENERGY_TS ), &    !-- cell center
        nProfile => size ( Profile, dim = 1 ) )

    allocate ( Radius ( nProfile + 1 ) )
    allocate ( RadialVelocity ( nProfile + 1 ) )
    allocate ( dRC ( nProfile ) )
    allocate ( RC ( nProfile ) )
    Radius ( 1 )          =  0.0_KDR
    RadialVelocity ( 1 )  =  0.0_KDR
    do iV = 2, nProfile + 1
      Radius         ( iV )  =  R ( iV - 1 )
      RadialVelocity ( iV )  =  V_R ( iV - 1 )
      dRC ( iV - 1 )  =  Radius ( iV )  -  Radius ( iV - 1 )
      RC  ( iV - 1 )  =  Radius ( iV - 1 )  +  0.5_KDR * dRC ( iV - 1 )
    end do

    allocate ( Density ( nProfile + 1 ) )
    allocate ( Temperature ( nProfile + 1 ) )
    allocate ( SpecificEnergy ( nProfile + 1 ) )
    allocate ( ElectronFraction ( nProfile + 1 ) )

    !-- First edge extrapolated
    Slope_N    =  ( N ( 2 )  -  N ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_T    =  ( T ( 2 )  -  T ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_SE   =  ( SE ( 2 )  -  SE ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_Y_E  =  ( Y_E ( 2 )  -  Y_E ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )

    Density          ( 1 )  =  N ( 1 )  &
                               +  Slope_N   * ( Radius ( 1 )  -  RC ( 1 ) )
    Temperature      ( 1 )  =  T ( 1 )  &
                               +  Slope_T   * ( Radius ( 1 )  -  RC ( 1 ) )
    SpecificEnergy   ( 1 )  =  SE ( 1 )  &
                               +  Slope_SE  * ( Radius ( 1 )  -  RC ( 1 ) )
    ElectronFraction ( 1 )  =  Y_E ( 1 )  &
                               +  Slope_Y_E * ( Radius ( 1 )  -  RC ( 1 ) )

    do iV = 2, nProfile + 1

      if ( iV <= nProfile ) then
        Slope_N    =  ( N ( iV )  -  N ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_T    =  ( T ( iV )  -  T ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_SE   =  ( SE ( iV )  -  SE ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_Y_E  =  ( Y_E ( iV )  -  Y_E ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
      else
        !-- Last edge extrapolated with same slope
      end if

      Density ( iV )  &
        =  N   ( iV - 1 )  +  Slope_N   * ( Radius ( iV )  -  RC ( iV - 1 ) )
      Temperature ( iV )  &
        =  T   ( iV - 1 )  +  Slope_T   * ( Radius ( iV )  -  RC ( iV - 1 ) )
      SpecificEnergy ( iV )  &
        =  SE  ( iV - 1 )  +  Slope_SE  * ( Radius ( iV )  -  RC ( iV - 1 ) )
      ElectronFraction ( iV )  &
        =  Y_E ( iV - 1 )  +  Slope_Y_E * ( Radius ( iV )  -  RC ( iV - 1 ) )

    end do !-- iV

    end associate !-- R, etc.

    call Show ( 'First few values' )
    call Show ( Profile ( 1 : 5, iRADIUS_TS ), 'RadiusTable' )
    call Show ( Radius ( 1 : 5 ), 'RadiusEdge' )
    call Show ( Profile ( 1 : 5, iRADIAL_VELOCITY_TS ), 'RadialVelocityTable' )
    call Show ( RadialVelocity ( 1 : 5 ), 'RadialVelocityEdge' )
    call Show ( Profile ( 1 : 5, iDENSITY_TS ), 'DensityTable' )
    call Show ( Density ( 1 : 5 ), 'DensityEdge' )
    call Show ( Profile ( 1 : 5, iTEMPERATURE_TS ), 'TemperatureTable' )
    call Show ( Temperature ( 1 : 5 ), 'TemperatureEdge' )
    call Show ( Profile ( 1 : 5, iSPECIFIC_ENERGY_TS ), 'SpecificEnergyTable' )
    call Show ( SpecificEnergy ( 1 : 5 ), 'SpecificEnergyEdge' )
    call Show ( Profile ( 1 : 5, iELECTRON_FRACTION_TS ), &
                'ElectronFractionTable' )
    call Show ( ElectronFraction ( 1 : 5 ), 'ElectronFractionEdge' )

    Radius          =  Radius          *  UNIT % CENTIMETER
    RadialVelocity  =  RadialVelocity  *  UNIT % CENTIMETER / UNIT % SECOND
    Density         =  Density         *  UNIT % MASS_DENSITY_CGS
    Temperature     =  Temperature     *  UNIT % KELVIN
    SpecificEnergy  =  SpecificEnergy  *  UNIT % ERG / UNIT % GRAM

    !-- SplineInterpolation initialization

    call SI ( iRADIAL_VELOCITY_SI ) % Initialize &
           ( Radius, RadialVelocity, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iDENSITY_SI ) % Initialize &
           ( Radius, Density, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iTEMPERATURE_SI ) % Initialize &
           ( Radius, Temperature, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iSPECIFIC_ENERGY_SI ) % Initialize &
           ( Radius, SpecificEnergy, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iELECTRON_FRACTION_SI ) % Initialize &
           ( Radius, ElectronFraction, VerbosityOption = CONSOLE % INFO_3 )

  end subroutine PrepareInterpolation


  subroutine ApplySourcesGravity ( S, Increment, Fluid, TimeStep )

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
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      KV_M_1, &
      KV_E, &
      D, &
      V_1, &
      dPhi_dR
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call ApplySourcesCurvilinear_Fluid_P ( S, Increment, Fluid, TimeStep )

    select type ( F => Fluid )
    class is ( Fluid_P_MHN_Form )

    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )

    select type ( Chart => S % Chart )
    class is ( Chart_SLD_Form )

    G => Chart % Geometry ( )

!     iStage = mod ( iStage, S % nStages ) + 1

!    if ( iStage == 1 ) &
      call ComputeGravity ( F, Chart % Atlas, G )

    call Chart % SetVariablePointer &
           ( Increment % Value ( :, iMomentum_1 ), KV_M_1 )
    call Chart % SetVariablePointer &
           ( Increment % Value ( :, iEnergy ), KV_E )
    call Chart % SetVariablePointer &
           ( F % Value ( :, F % CONSERVED_DENSITY ), D )
    call Chart % SetVariablePointer &
           ( F % Value ( :, F % VELOCITY_U ( 1 ) ), V_1 )
    call Chart % SetVariablePointer &
           ( Acceleration_C, dPhi_dR )

    call ApplySourcesKernel &
           ( KV_M_1, KV_E, D, V_1, dPhi_dR, CONSTANT % GRAVITATIONAL, &
             TimeStep, oV = Chart % nGhostLayers )

    end select !-- Chart
    end select !-- F

    nullify ( KV_M_1, KV_E, dPhi_dR, D, V_1 )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplySourcesGravity


  subroutine ComputeGravity ( F, PS, G )

    class ( CurrentTemplate ), intent ( in ) :: &
      F
    class ( AtlasHeaderForm ), intent ( in ) :: &
      PS
    type ( GeometryFlatForm ), intent ( in ) :: &
      G

    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Phi_C, &
      dPhi_dR_C, &
      R

    select type ( F )
    class is ( Fluid_P_MHN_Form )

    select type ( PS )
    class is ( Atlas_SC_Form )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )

    call ComputeGravityEdge &
           ( C, F % Value ( :, F % CONSERVED_DENSITY ), &
             G % Value ( :, G % VOLUME_JACOBIAN ), &
             G % VALUE ( :, G % CENTER ( 1 ) ), &
             G % VALUE ( :, G % WIDTH ( 1 ) ), CONSTANT % GRAVITATIONAL, &
             Potential )

    call Clear ( Potential_C )
    call Clear ( Acceleration_C )
    call C % SetVariablePointer ( Potential_C,    Phi_C )
    call C % SetVariablePointer ( Acceleration_C, dPhi_dR_C )
    call C % SetVariablePointer ( G % VALUE ( :, G % CENTER ( 1 ) ), R )

    call ComputeGravityCenter ( Phi_C, dPhi_dR_C, C, R, Potential )

    ! !-- 1D

    ! call Clear ( Potential_C )

    ! oV   =  Chart % nGhostLayers ( 1 )
    ! nV   =  Chart % nCellsBrick ( 1 )
    ! oVM  =  ( Chart % iaBrick ( 1 ) - 1 ) * Chart % nCellsBrick ( 1 ) &
    !         -  Chart % nGhostLayers ( 1 )
    ! do iV = oV + 1, oV + nV
    !   Potential_C ( iV )  &
    !     =  0.5_KDR * ( Potential ( oVM + iV - 1 )  +  Potential ( oVM + iV ) )
    ! end do
    ! Potential_C ( oV + nV + 1 )  =  Potential ( oVM + oV + nV )

    end select !-- C
    end select !-- PS
    end select !-- F

    nullify ( Phi_C, dPhi_dR_C, R )

  end subroutine ComputeGravity


  subroutine ApplySourcesKernel ( KV_M_1, KV_E, D, V_1, dPhi_dR, G, dT, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      KV_M_1, &
      KV_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      D, &
      V_1, &
      dPhi_dR
    real ( KDR ), intent ( in ) :: &
      G, &
      dT
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      oV

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      lV, uV
    real ( KDR ) :: &
      dPhi_dR_C

    lV = 1
    where ( shape ( D ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( D ) > 1 )
      uV = shape ( D ) - oV
    end where

    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
          KV_M_1 ( iV, jV, kV )  &
            =  KV_M_1 ( iV, jV, kV )  &
               -  dT  *  D ( iV, jV, kV )  *  dPhi_dR ( iV, jV, kV )
          KV_E ( iV, jV, kV )  &
            =  KV_E ( iV, jV, kV )  &
               -  dT  *  D ( iV, jV, kV )  *  dPhi_dR ( iV, jV, kV ) &
                      *  V_1 ( iV, jV, kV ) 
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ApplySourcesKernel


  subroutine ComputeGravityEdge ( C, D, VJ, R, dR, G, Phi )

    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      D, &
      VJ, &
      R, &
      dR
    real ( KDR ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( 0 : ), intent ( out ) :: &
      Phi

    integer ( KDI ) :: &
      iV, &
      lV, &
      uV
    real ( KDR ), dimension ( 0 : size ( Phi ) - 1 ) :: &
      SHR, &  !-- SolidHarmonicRegular
      SHI, &  !-- SolidHarmonicIrregular
      R_I     !-- R_Inner
    real ( KDR ), dimension ( size ( Phi ( 1 : ) ) ) :: &
      D_Ave, &
      dV, &
      R_C  !-- R_Center
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select case ( C % nDimensions )
    case ( 1 )

      call CO % Initialize &
             ( C % Atlas % Communicator, &
               nOutgoing = [ C % nCellsBrick ( 1 ) ], &
               nIncoming = [ C % nCells ( 1 ) ] )

      CO % Outgoing % Value = pack ( D, mask = C % IsProperCell )
      call CO % Gather ( )
      D_Ave  =  CO % Incoming % Value

      CO % Outgoing % Value = pack ( VJ * dR, mask = C % IsProperCell )
      call CO % Gather ( )
      dV  =  CO % Incoming % Value

      CO % Outgoing % Value = pack ( R, mask = C % IsProperCell )
      call CO % Gather ( )
      R_C  =  CO % Incoming % Value

      CO % Outgoing % Value = pack ( dR, mask = C % IsProperCell )
      call CO % Gather ( )
      lV = 0
      uV = size ( Phi ) - 1
      R_I ( lV : uV - 1 )   =  R_C  -  0.5_KDR  *  CO % Incoming % Value
      R_I ( uV ) =  R_C ( size ( R_C ) )  &
                    +  0.5_KDR  *  CO % Incoming % Value ( size ( R_C ) )

    case default
      call Show ( 'Dimensionality not implemented', CONSOLE % ERROR )
      call Show ( 'WoosleyHeger_07_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeGravityalPotentialKernel', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- nDimensions

    SHR ( 0 )  =  0.0_KDR
    do iV = 1, size ( SHR ( 1 : ) )
      SHR ( iV )  =  SHR ( iV - 1 )  +  D_Ave ( iV ) * dV ( iV )
    end do !-- iV

    SHI ( size ( SHI ) - 1 )  =  0.0_KDR
!call Show ( size ( SHI ) - 1, '>>> size ( SHI ) - 1' )
    do iV = size ( SHI ) - 2, 0, -1
!call Show ( iV, '>>> iV' )
      SHI ( iV )  =  SHI ( iV + 1 )  &
                     +  D_Ave ( iV + 1 ) / R_C ( iV + 1 ) * dV ( iV + 1 )  
    end do !-- iV
!call Show ( SHI, 'SHI' )

    Phi ( 0 )  =  - SHI ( 0 )
    do iV  =  1, size ( Phi ) - 1
      Phi ( iV )  =  - G * ( SHR ( iV ) / R_I ( iV )  +  SHI ( iV ) )
    end do !-- iV

  end subroutine ComputeGravityEdge


  subroutine ComputeGravityCenter ( Phi_C, dPhi_dR_C, C, R, Phi )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      Phi_C, &
      dPhi_dR_C
    class ( Chart_SLD_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      R
    real ( KDR ), dimension ( 0 : ), intent ( in ) :: &
      Phi

    integer ( KDI ) :: &
      iV, jV, kV, &
      oV, &
      oVM
    integer ( KDI ), dimension ( 3 ) :: &
      lV, uV

    lV = 1
    where ( shape ( Phi_C ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( Phi_C ) > 1 )
      uV = shape ( Phi_C ) - oV
    end where

    oV   =  C % nGhostLayers ( 1 )
    oVM  =  ( C % iaBrick ( 1 ) - 1 ) * C % nCellsBrick ( 1 ) &
            -  C % nGhostLayers ( 1 )

    !-- Cell-centered Phi
    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
          Phi_C ( iV, jV, kV )  &
            =  0.5_KDR * ( Phi ( oVM + iV - 1 )  +  Phi ( oVM + iV ) )
        end do
      end do
    end do
    !$OMP end parallel do

    !-- Radial boundary and ghost values of cell-centered Phi 
    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )

        iV = lV ( 1 ) - 1
        if ( C % iaBrick ( 1 ) == 1 ) then 
          !-- inner boundary: vanishing gradient
          Phi_C ( iV, jV, kV )  =  Phi_C ( iV + 1, jV, kV )
        else
          !-- ghost cell, interior value
          Phi_C ( iV, jV, kV ) &
            =  0.5_KDR * ( Phi ( oVM + iV - 1 )  +  Phi ( oVM + iV ) )
        end if
       
        iV = uV ( 1 ) + 1
        if ( C % iaBrick ( 1 ) == C % nBricks ( 1 ) ) then
          !-- outer boundary: M / r
          Phi_C ( iV, jV, kV ) &
            =  Phi_C ( iV - 1, jV, kV )  *  R ( iV - 1, jV, kV ) &
                                         /  R ( iV, jV, kV )
        else
          !-- ghost cell, interior value
          Phi_C ( iV, jV, kV ) &
            =  0.5_KDR * ( Phi ( oVM + iV - 1 )  +  Phi ( oVM + iV ) )
        end if

      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          dPhi_dR_C ( iV, jV, kV ) &
            =  ( Phi_C ( iV + 1, jV, kV )  -  Phi_C ( iV - 1, jV, kV ) ) &
               / ( R ( iV + 1, jV, kV )  -  R ( iV - 1, jV, kV ) )

        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeGravityCenter


!   function LocalMax ( IsProperCell, V ) result ( ML ) 

!     logical ( KDL ), dimension ( : ), intent ( in ) :: &
!       IsProperCell
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       V
!     real ( KDR ) :: &
!       ML

!     integer ( KDI ) :: &
!       iV

!     ML = - huge ( 0.0_KDR )
!     !$OMP parallel do private ( iV ) &
!     !$OMP reduction ( max : ML )
!     do iV = 1, size ( V )
!       if ( IsProperCell ( iV ) ) &
!         ML  =  max ( ML, V ( iV ) )
!     end do !-- iV
!     !$OMP end parallel do
 
!   end function LocalMax


end module WoosleyHeger_07_Header_Form
