module WoosleyHeger_07__Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  use Fluid_P_MHN__Form
  use Fluid_ASC__Form
  use Tally_F_P__Form
  use Diagnostics_WH07__Form
  use Diagnostics_WH07_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: WoosleyHeger_07_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &  !-- 3
      SetWriteTimeInterval
    final :: &
      Finalize  
  end type WoosleyHeger_07_Form

    private :: &
      SetFluid, &
      ApplySources

      private :: &
        PrepareInterpolation, &
        ApplySourcesKernel, &
        LocalMax

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

    integer ( KDI ), private :: &
      iStage = 0
    real ( KDR ), dimension ( : ), allocatable :: &
      EnclosedMass, &
      Potential
    type ( Diagnostics_WH07_ASC_Form ), dimension ( : ), allocatable :: &
      Diagnostics_ASC

contains


  subroutine Initialize ( WH, Name )

    class ( WoosleyHeger_07_Form ), intent ( inout ) :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iS,  & !-- iStage
      nCellsCore, &
      nCellsRadius, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax, &
      CellRatio, &
      FinishTime, &
      CourantFactor
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
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
      CoordinateUnit, &
      VelocityUnit
    character ( 1 + 2 ) :: &
      StageNumber
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: WH % PositionSpace )
    select type ( PS => WH % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    CoordinateSystem = 'SPHERICAL'
    CoordinateUnit   = [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

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

    ! CellRatio = 1.026_KDR  !-- nCells ( 1 ) = 256
    ! call PROGRAM_HEADER % GetParameter ( CellRatio, 'CellRatio' )

    ! Spacing = [ 'GEOMETRIC', 'EQUAL    ', 'EQUAL    ' ]
    ! Ratio   = [ CellRatio, 0.0_KDR, 0.0_KDR ]

    ! call PS % CreateChart &
    !        ( SpacingOption = Spacing, &
    !          CoordinateSystemOption = CoordinateSystem, &
    !          CoordinateUnitOption = CoordinateUnit, &
    !          MinCoordinateOption = MinCoordinate, &
    !          MaxCoordinateOption = MaxCoordinate, &
    !          RatioOption = Ratio, &
    !          nCellsOption = nCells )

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    RadiusCore = 40.0_KDR  *  UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( RadiusCore, 'RadiusCore' )

    nCellsCore = 128  !-- Number of central cells with equal spacing
    call PROGRAM_HEADER % GetParameter ( nCellsCore, 'nCellsCore' )

    nCellsRadius = 6.5 * nCellsCore  !-- Aiming for roughly 10,000 km
    call PROGRAM_HEADER % GetParameter ( nCellsRadius, 'nCellsRadius' )

    nCellsPolar     = 3 * nCellsCore
    nCellsAzimuthal = 2 * nCellsPolar

    call Show ( 'Mesh core parameters' )
    call Show ( RadiusCore, UNIT % KILOMETER, 'RadiusCore' )
    call Show ( nCellsCore, 'nCellsCore' )
    call Show ( RadiusCore / nCellsCore, UNIT % KILOMETER, 'CellWidthCore' )

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
             CoordinateUnitOption = CoordinateUnit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nEqualOption = nCellsCore )

    !-- Geometry of PositionSpace

    allocate ( WH % Geometry_ASC )
    associate ( GA => WH % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Fluid

    TimeUnit = UNIT % SECOND

    VelocityUnit       =  CoordinateUnit / TimeUnit 
    MassDensityUnit    =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit  =  UNIT % MASS_DENSITY_CGS  &
                            *  UNIT % SPEED_OF_LIGHT ** 2
    TemperatureUnit    =  UNIT % MEV

    MassUnit             =  UNIT % SOLAR_MASS
    EnergyUnit           =  MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    MomentumUnit         =  MassUnit  *  UNIT % SPEED_OF_LIGHT
    AngularMomentumUnit  =  CoordinateUnit ( 1 ) &
                            *  MassUnit  *  UNIT % SPEED_OF_LIGHT
    
    allocate ( Fluid_ASC_Form :: WH % Current_ASC )
    select type ( FA => WH % Current_ASC )
    class is ( Fluid_ASC_Form )

    call FA % Initialize &
           ( PS, 'MEAN_HEAVY_NUCLEUS', &
             VelocityUnitOption = VelocityUnit, &
             MassDensityUnitOption = MassDensityUnit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             TemperatureUnitOption   = TemperatureUnit, &
             MassUnitOption = MassUnit, EnergyUnitOption = EnergyUnit, &
             MomentumUnitOption = MomentumUnit, &
             AngularMomentumUnitOption = AngularMomentumUnit )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: WH % Step )
    select type ( S => WH % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources

    !-- Diagnostics

    allocate ( Diagnostics_ASC ( S % nStages ) )
    do iS = 1, S % nStages
      write ( StageNumber, fmt = '(a1,i2.2)' ) '_', iS
      associate ( DA => Diagnostics_ASC ( iS ) )
      call DA % Initialize ( PS, 'Diagnostics' // StageNumber )
      end associate !-- DA
    end do !-- iS

    !-- Set fluid and initialize Integrator template

    FinishTime    = 1.0_KDR * UNIT % SECOND
    CourantFactor = 0.7_KDR

    call SetFluid ( WH )
    call WH % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, FinishTimeOption = FinishTime, &
             CourantFactorOption = CourantFactor, nWriteOption = 30 )

    !-- Cleanup

    end select !-- S
    end select !-- FA
    end select !-- PS

  end subroutine Initialize


  subroutine SetWriteTimeInterval ( I )

    class ( WoosleyHeger_07_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iProcess, &
      iRadius
    real ( KDR ) :: &
      VelocityMax, &
      VelocityMaxRadius, &
      DensityAve, &
      TimeScaleDensityAve, &
      TimeScaleVelocityMax
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_MHN_Form ), pointer :: &
      F

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_MHN ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( Grid => PS % Chart )
    class is ( Chart_SL_Template )

    associate ( C => PS % Communicator ) 

    !-- Find max velocity
    allocate ( CO )
    call CO % Initialize ( C, nOutgoing = [ 1 ], nIncoming = [ C % Size ] )
    CO % Outgoing % Value ( 1 ) &
      =  LocalMax ( Grid % IsProperCell, &
                    abs ( F % Value ( :, F % VELOCITY_U ( 1 ) ) ) ) 
    call CO % Gather ( )
    VelocityMax  =  maxval ( CO % Incoming % Value )
    iProcess     =  maxloc ( CO % Incoming % Value, dim = 1 )  -  1
    deallocate ( CO )

    !-- Find radius of max velocity
    allocate ( CO )
    call CO % Initialize &
           ( C, nOutgoing = [ 1 ], nIncoming = [ 1 ], RootOption = iProcess )
    if ( C % Rank == iProcess ) then
      iRadius  =  maxloc ( abs ( F % Value ( :, F % VELOCITY_U ( 1 ) ) ), &
                           dim = 1 )
      CO % Outgoing % Value ( 1 )  =  G % Value ( iRadius, G % CENTER ( 1 ) )
    end if
    call CO % Broadcast ( )
    VelocityMaxRadius = CO % Incoming % Value ( 1 )
    deallocate ( CO )

    select type ( TI => FA % TallyInterior )
    class is ( Tally_F_P_Form )
    DensityAve  =  TI % Value ( TI % BARYON_NUMBER ) / VelocityMaxRadius ** 3 
    end select !-- TI

    TimeScaleVelocityMax &
      =  VelocityMaxRadius  /  VelocityMax
    TimeScaleDensityAve &
      =  ( CONSTANT % GRAVITATIONAL  *  DensityAve ) ** ( -0.5_KDR )

    I % WriteTimeInterval  &
      =  min ( TimeScaleDensityAve, TimeScaleVelocityMax )  /  I % nWrite

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( VelocityMax, Grid % CoordinateUnit ( 1 ) / I % TimeUnit, &
                'VelocityMax', I % IGNORABILITY )
    call Show ( VelocityMaxRadius, Grid % CoordinateUnit ( 1 ), &
                'VelocityMaxRadius', I % IGNORABILITY )
    call Show ( DensityAve, FA % MassDensityUnit, 'DensityAve', &
                I % IGNORABILITY )
    call Show ( TimeScaleDensityAve, I % TimeUnit, 'TimeScaleDensityAve', &
                I % IGNORABILITY )
    call Show ( TimeScaleVelocityMax, I % TimeUnit, 'TimeScaleVelocityMax', &
                I % IGNORABILITY )

    end associate !-- Comm
    end select !-- C
    end select !-- PS
    end select !-- FA
    nullify ( G, F )

  end subroutine SetWriteTimeInterval


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_Form ), intent ( inout ) :: &
      WH

    call WH % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetFluid ( WH )

    type ( WoosleyHeger_07_Form ), intent ( inout ) :: &
      WH

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      SE  !-- SpecificEnergy
    type ( SplineInterpolationForm ), dimension ( 5 ) :: &
      SI
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_MHN_Form ), pointer :: &
      F

    call PrepareInterpolation ( SI )

    select type ( FA => WH % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_MHN ( )

    select type ( PS => WH % PositionSpace )
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

    end do

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    call F % ComputeFromPrimitive ( G )

    end associate !-- N, etc.
    end select !-- PS
    end select !-- FA
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
!iV, &
      iMomentum_1, &
      iEnergy
! real ( KDR ), dimension ( Increment % nValues ) :: &
!   Force, &
!   Power
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      KV_M_1, &
      KV_E, &
      F_G_M, &
      F_G_Phi, &
      F_G_Phi_VJ, &
      Phi_C, &
      F_G_Phi_C, &
      D, &
      V_1, &
      R, &
      dR, &
      VJ, &
      VJ_I, &
      dLVJ_dX1
    type ( VariableGroupForm ) :: &
      G_I  !-- GeometryInner
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Diagnostics_WH07_Form ), pointer :: &
      D_WH
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

    iStage = mod ( iStage, S % nStages ) + 1

    if ( .not. allocated ( EnclosedMass ) ) then
      allocate ( EnclosedMass ( 0 : Chart % nCells ( 1 ) ) )  !-- edges
      allocate ( Potential ( 0 : Chart % nCells ( 1 ) ) )  !-- edges
    end if

!    if ( iStage == 1 ) &
      call ComputeGravitationalPotential &
             ( Chart, F % Value ( :, F % CONSERVED_DENSITY ), &
               G % Value ( :, G % VOLUME_JACOBIAN ), &
               G % VALUE ( :, G % CENTER ( 1 ) ), &
               G % VALUE ( :, G % WIDTH ( 1 ) ), CONSTANT % GRAVITATIONAL, &
               EnclosedMass, Potential )

!    call Show ( M, UNIT % SOLAR_MASS, '>>> M' )
!    call Show ( Phi, UNIT % SPEED_OF_LIGHT ** (-2), '>>> Phi' )
    
    call G_I % Initialize ( shape ( G % Value ) )
    call G % ComputeReconstruction ( G_I, Chart % nDimensions, iDimension = 1 )

    D_WH => Diagnostics_ASC ( iStage ) % Diagnostics_WH07 ( )
    associate &
      ( DV_F_P => D_WH % Value ( :, D_WH % PRESSURE_FORCE ), &
        DV_F_N => D_WH % Value ( :, D_WH % NET_FORCE ), &
        DV_P_N => D_WH % Value ( :, D_WH % NET_POWER ) )

    DV_F_P = Increment % Value ( :, iMomentum_1 ) / TimeStep

    call Chart % SetVariablePointer &
           ( Increment % Value ( :, iMomentum_1 ), KV_M_1 )
    call Chart % SetVariablePointer &
           ( Increment % Value ( :, iEnergy ), KV_E )
    call Chart % SetVariablePointer &
           ( D_WH % Value ( :, D_WH % GRAVITATIONAL_FORCE_M ), F_G_M )
    call Chart % SetVariablePointer &
           ( D_WH % Value ( :, D_WH % GRAVITATIONAL_FORCE_PHI ), F_G_Phi )
    call Chart % SetVariablePointer &
           ( D_WH % Value ( :, D_WH % GRAVITATIONAL_FORCE_PHI_VJ ), F_G_Phi_VJ)
    call Chart % SetVariablePointer &
           ( D_WH % Value ( :, D_WH % GRAVITATIONAL_POTENTIAL_C ), Phi_C )
    call Chart % SetVariablePointer &
           ( D_WH % Value ( :, D_WH % GRAVITATIONAL_FORCE_PHI_C ), F_G_Phi_C )
    call Chart % SetVariablePointer &
           ( F % Value ( :, F % CONSERVED_DENSITY ), D )
    call Chart % SetVariablePointer &
           ( F % Value ( :, F % VELOCITY_U ( 1 ) ), V_1 )
    call Chart % SetVariablePointer &
           ( G % Value ( :, G % CENTER ( 1 ) ), R )
    call Chart % SetVariablePointer &
           ( G % Value ( :, G % WIDTH ( 1 ) ), dR )
    call Chart % SetVariablePointer &
           ( G % Value ( :, G % VOLUME_JACOBIAN ), VJ )
    call Chart % SetVariablePointer &
           ( G_I % Value ( :, G % VOLUME_JACOBIAN ), VJ_I )
    call Chart % SetVariablePointer &
           ( S % dLogVolumeJacobian_dX ( 1 ) % Value, dLVJ_dX1 )
    call ApplySourcesKernel &
           ( KV_M_1, KV_E, F_G_M, F_G_Phi, F_G_Phi_VJ, Phi_C, F_G_Phi_C, &
             Chart, EnclosedMass, Potential, D, V_1, R, dR, VJ, VJ_I, &
             dLVJ_dX1, CONSTANT % GRAVITATIONAL, TimeStep, &
             oV = Chart % nGhostLayers, &
             oVM = ( Chart % iaBrick ( 1 ) - 1 ) * Chart % nCellsBrick ( 1 ) &
                   -  Chart % nGhostLayers ( 1 ) )
    ! call ApplySourcesKernel &
    !        ( Increment % Value ( :, iMomentum_1 ), &
    !          Increment % Value ( :, iEnergy ), &
    !          Chart % IsProperCell, M, &
    !          F % Value ( :, F % CONSERVED_DENSITY ), &
    !          F % Value ( :, F % VELOCITY_U ( 1 ) ), &
    !          G % Value ( :, G % CENTER ( 1 ) ), &
    !          CONSTANT % GRAVITATIONAL, TimeStep ) 

! call Chart % ExchangeGhostData ( Increment )
! Force  =  Increment % Value ( :, iMomentum_1 )
! Power  =  Increment % Value ( :, iEnergy )
! do iV = 2, Increment % nValues - 1
!   if ( PT ( iV ) > 0.0_KDR .and. iV > 20 ) then
!     Increment % Value ( iV, iMomentum_1 ) &
!       =  ( Force ( iV - 1 ) + Force ( iV + 1 ) ) / 2.0_KDR
!     Increment % Value ( iV, iEnergy ) &
!       =  ( Power ( iV - 1 ) + Power ( iV + 1 ) ) / 2.0_KDR
!   end if
! end do !-- iV

    DV_F_N = Increment % Value ( :, iMomentum_1 ) / TimeStep
    DV_P_N = Increment % Value ( :, iEnergy ) / TimeStep

    end associate !-- DV_F_P, etc.
    end select !-- Chart
    end select !-- F

    nullify ( KV_M_1, KV_E, F_G_M, F_G_Phi, F_G_Phi_VJ, Phi_C, F_G_Phi_C, &
              D, V_1, R, G, D_WH )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplySources

  
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

    end do

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


  subroutine ComputeGravitationalPotential &
               ( C, D, VJ, R, dR, G, M, Phi )

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
      M, &
      Phi

    integer ( KDI ) :: &
      iV, &
      lV, &
      uV
    real ( KDR ), dimension ( 0 : size ( M ) - 1 ) :: &
      SHI, &  !-- SolidHarmonicIrregular
              !-- Enclosed mass M is SolidHarmonicRegular
      R_I     !-- R_Inner
    real ( KDR ), dimension ( size ( M ( 1 : ) ) ) :: &
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
      uV = size ( M ) - 1
      R_I ( lV : uV - 1 )   =  R_C  -  0.5_KDR  *  CO % Incoming % Value
      R_I ( uV ) =  R_C ( size ( R_C ) )  &
                    +  0.5_KDR  *  CO % Incoming % Value ( size ( R_C ) )

    case default
      call Show ( 'Dimensionality not implemented', CONSOLE % ERROR )
      call Show ( 'WoosleyHeger_07_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeGravitationalPotential', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- nDimensions

!    call Show ( N, UNIT % MASS_DENSITY_CGS, '>>> N' )
!    call Show ( N_Ave, UNIT % MASS_DENSITY_CGS, '>>> N_Ave' )

    M ( 0 )  =  0.0_KDR
    do iV = 1, size ( M ( 1 : ) )
      M ( iV )  =  M ( iV - 1 )  +  D_Ave ( iV ) * dV ( iV )
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
      Phi ( iV )  =  - G * ( M ( iV ) / R_I ( iV )  +  SHI ( iV ) )
    end do !-- iV

  end subroutine ComputeGravitationalPotential


  subroutine ApplySourcesKernel &
               ( KV_M_1, KV_E, F_G_M, F_G_Phi, F_G_Phi_VJ, Phi_C, F_G_Phi_C, &
                 Chart, M, Phi, D, V_1, R, dR, VJ, VJ_I, dLVJ_dX1, &
                 G, dT, oV, oVM )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      KV_M_1, &
      KV_E, &   
      F_G_M, &
      F_G_Phi, &
      F_G_Phi_VJ, &
      Phi_C, &
      F_G_Phi_C
    class ( Chart_SLD_Form ), intent ( in ) :: &
      Chart
    real ( KDR ), dimension ( 0 : ), intent ( in ) :: &
      M, &
      Phi
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      D, &
      V_1, &
      R, &
      dR, &
      VJ, &
      VJ_I, &
      dLVJ_dX1
    real ( KDR ) :: &
      G, &
      dT
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      oV
    integer ( KDI ), intent ( in ) :: &
      oVM

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      lV, uV, &
      nV
    real ( KDR ) :: &
      M_C, &  !-- M_Center
      dPhi, &
      VJ_Ave, &
      dPhi_dR_C, &
      DivPhi, &
      dPhi_dR_Div

    lV = 1
    where ( shape ( D ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( D ) > 1 )
      uV = shape ( D ) - oV
    end where

    nV = shape ( D )
    
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
        if ( Chart % iaBrick ( 1 ) == 1 ) then 
          !-- inner boundary: vanishing gradient
          Phi_C ( iV, jV, kV )  =  Phi_C ( iV + 1, jV, kV )
        else
          !-- ghost cell, interior value
          Phi_C ( iV, jV, kV ) &
            =  0.5_KDR * ( Phi ( oVM + iV - 1 )  +  Phi ( oVM + iV ) )
        end if
       
        iV = uV ( 1 ) + 1
        if ( Chart % iaBrick ( 1 ) == Chart % nBricks ( 1 ) ) then
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
! call Show ( oVM, '>>> oVM' )
! call Show ( iV, '>>> iV' )
! call Show ( oVM + iV - 1, '>>> oVM + iV - 1' )
! call Show ( oVM + iV, '>>> oVM + iV' )

          M_C  =  0.5_KDR * ( M ( oVM + iV - 1 )  +  M ( oVM + iV ) )
          F_G_M ( iV, jV, kV )  =  - G * M_C * D ( iV, jV, kV )  &
                                     /  R ( iV, jV, kV ) ** 2

          ! dPhi  =  Phi ( oVM + iV )  -  Phi ( oVM + iV - 1 )
          ! F_G_Phi ( iV, jV, kV )  =  - D ( iV, jV, kV )  *  dPhi  &
          !                              /  dR ( iV, jV, kV )

          ! VJ_Ave  =  0.5_KDR * ( VJ_I ( iV, jV, kV ) &
          !                        +  VJ_I ( iV + 1, jV, kV ) )
          ! F_G_Phi_VJ ( iV, jV, kV )  &
          !   =  F_G_Phi ( iV, jV, kV )  *  VJ_Ave  /  VJ ( iV, jV, kV )

          DivPhi  =  (    VJ_I ( iV + 1, jV, kV ) * Phi ( oVM + iV ) &
                       -  VJ_I ( iV, jV, kV ) * Phi ( oVM + iV - 1 ) ) &
                     / ( VJ ( iV, jV, kV )  *  dR ( iV, jV, kV ) )

          dPhi_dR_Div  &
            =  DivPhi  -  dLVJ_dX1 ( iV, jV, kV ) * Phi_C ( iV, jV, kV )

          F_G_Phi_VJ ( iV, jV, kV )  =  - D ( iV, jV, kV )  *  dPhi_dR_Div

          dPhi_dR_C  &
            =  ( Phi_C ( iV + 1, jV, kV )  -  Phi_C ( iV - 1, jV, kV ) ) &
               / ( R ( iV + 1, jV, kV )  -  R ( iV - 1, jV, kV ) )

          F_G_Phi_C ( iV, jV, kV )  =  - D ( iV, jV, kV )  *  dPhi_dR_C

!call Show ( iV, '>>> iV' )
!call show ( [ F_G_M ( iV, jV, kV ), F_G_Phi ( iV, jV, kV ), &
!              F_G_Phi_VJ ( iV, jV, kV ) ], &
!            '>>> F_G_M, F_G_Phi, F_G_Phi_VJ' )
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
          KV_M_1 ( iV, jV, kV )  &
            =  KV_M_1 ( iV, jV, kV )  +  dT * F_G_Phi_C ( iV, jV, kV )
          KV_E ( iV, jV, kV )  &
            =  KV_E ( iV, jV, kV )  &
               +  dT * F_G_Phi_C ( iV, jV, kV ) * V_1 ( iV, jV, kV ) 
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ApplySourcesKernel


  function LocalMax ( IsProperCell, V ) result ( ML ) 

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V
    real ( KDR ) :: &
      ML

    integer ( KDI ) :: &
      iV

    ML = - huge ( 0.0_KDR )
    !$OMP parallel do private ( iV ) &
    !$OMP reduction ( max : ML )
    do iV = 1, size ( V )
      if ( IsProperCell ( iV ) ) &
        ML  =  max ( ML, V ( iV ) )
    end do !-- iV
    !$OMP end parallel do
 
  end function LocalMax


end module WoosleyHeger_07__Form
