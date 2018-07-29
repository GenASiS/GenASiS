module FluidCentral_Template

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ApplyGravity_F__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    FluidCentralTemplate
      type ( MeasuredValueForm ), dimension ( 3 ) :: &
        CoordinateUnit
      logical ( KDL ) :: &
        Dimensionless
      type ( GridImageStreamForm ), allocatable :: &
        GridImageStream_SA_EB
      type ( Atlas_SC_Form ), allocatable :: &
        PositionSpace_SA, &  !-- SphericalAverage
        PositionSpace_EB     !-- EnclosedBaryons
      type ( Fluid_ASC_Form ), allocatable :: &
        Fluid_ASC_SA, &  !-- SphericalAverage
        Fluid_ASC_EB     !-- EnclosedBaryons
  contains
    procedure, public, pass :: &
      InitializeTemplate_FC
    procedure ( IPS ), private, pass, deferred :: &
      InitializePositionSpace
    procedure ( IG ), private, pass, deferred :: &
      InitializeGeometry
    procedure, public, pass :: &  !-- 2
      OpenGridImageStreams
    procedure, public, pass :: &  !-- 2
      OpenManifoldStreams
    procedure, public, pass :: &  !-- 3
      Write
    procedure, public, pass :: &
      FinalizeTemplate_FC
  end type FluidCentralTemplate

  abstract interface

    subroutine IPS ( FC, RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
                     RadialRatioOption, nCellsCoreOption, nCellsPolarOption )
      use Basics
      import FluidCentralTemplate
      class ( FluidCentralTemplate ), intent ( inout ) :: &
        FC
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusMinOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsCoreOption, &
        nCellsPolarOption
    end subroutine IPS

    subroutine IG ( FC, GA, PS, GeometryType, CentralMassOption )
      use Basics
      use Mathematics
      use Spaces
      import FluidCentralTemplate
      class ( FluidCentralTemplate ), intent ( inout ) :: &
        FC
      type ( Geometry_ASC_Form ), intent ( inout ) :: &
        GA
      class ( Atlas_SC_Form ), intent ( in ) :: &
        PS
      character ( * ), intent ( in )  :: &
        GeometryType
      real ( KDR ), intent ( in ), optional :: &
        CentralMassOption
    end subroutine IG

  end interface

    private :: &
      InitializeDiagnostics, &
      ComputeSphericalAverage, &
      ComputeEnclosedBaryons, &
      ApplySources

contains


  subroutine InitializeTemplate_FC &
               ( FC, Name, FluidType, GeometryType, DimensionlessOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, &
                 RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
                 RadialRatioOption, CentralMassOption, nWriteOption, &
                 nCellsCoreOption, nCellsPolarOption )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      LimiterParameterOption, &
      ShockThresholdOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusMinOption, &
      RadialRatioOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption, &
      nCellsCoreOption, &
      nCellsPolarOption

    real ( KDR ) :: &
      FinishTime
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_D_Unit
    class ( Fluid_D_Form ), pointer :: &
      F

    if ( FC % Type == '' ) &
      FC % Type = 'a FluidCentral'

    FC % Dimensionless = .false.
    if ( present ( DimensionlessOption ) ) &
      FC % Dimensionless = DimensionlessOption

    if ( .not. FC % Dimensionless ) then
      FC % CoordinateUnit  &
        =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]
      TimeUnit &
        = UNIT % SECOND
    end if

    if ( present ( TimeUnitOption ) ) &
      TimeUnit = TimeUnitOption


    !-- PositionSpace

    call FC % InitializePositionSpace &
           ( RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
             RadialRatioOption, nCellsCoreOption, nCellsPolarOption )

    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    call FC % InitializeGeometry ( GA, PS, GeometryType, CentralMassOption )

    call PS % SetGeometry ( GA )


    !-- Fluid

    allocate ( Fluid_ASC_Form :: FC % Current_ASC )
    select type ( FA => FC % Current_ASC )
    class is ( Fluid_ASC_Form )

    if ( FC % Dimensionless ) then
      call FA % Initialize ( PS, FluidType )
    else
      BaryonMassUnit                =  UNIT % ATOMIC_MASS_UNIT
      NumberDensityUnit             =  UNIT % NUMBER_DENSITY_NUCLEAR
      EnergyDensityUnit             =  UNIT % ENERGY_DENSITY_NUCLEAR
      TemperatureUnit               =  UNIT % MEGA_ELECTRON_VOLT
      Velocity_U_Unit               =  FC % CoordinateUnit  /  TimeUnit
      MomentumDensity_D_Unit        =  BaryonMassUnit * NumberDensityUnit &
                                       * Velocity_U_Unit
      MomentumDensity_D_Unit ( 2 )  =  MomentumDensity_D_Unit ( 2 ) &
                                       *  FC % CoordinateUnit ( 1 ) ** 2
      MomentumDensity_D_Unit ( 3 )  =  MomentumDensity_D_Unit ( 3 ) &
                                       *  FC % CoordinateUnit ( 1 ) ** 2
      call FA % Initialize &
             ( PS, FluidType, &
               Velocity_U_UnitOption         =  Velocity_U_Unit, &
               MomentumDensity_D_UnitOption  =  MomentumDensity_D_Unit, &
               BaryonMassUnitOption          =  BaryonMassUnit, &
               NumberDensityUnitOption       =  NumberDensityUnit , &
               EnergyDensityUnitOption       =  EnergyDensityUnit, &
               TemperatureUnitOption         =  TemperatureUnit, &
               NumberUnitOption              =  UNIT % SOLAR_BARYON_NUMBER, &
               EnergyUnitOption              =  UNIT % ENERGY_SOLAR_MASS, &
               MomentumUnitOption            =  UNIT % MOMENTUM_SOLAR_MASS, &
               AngularMomentumUnitOption     =  UNIT % SOLAR_KERR_PARAMETER, &
               TimeUnitOption                =  TimeUnit, &
               BaryonMassReferenceOption     =  CONSTANT % ATOMIC_MASS_UNIT, &
               LimiterParameterOption        =  LimiterParameterOption, &
               ShockThresholdOption          =  ShockThresholdOption )
    end if

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FC % Step )
    select type ( S => FC % Step )
    class is ( Step_RK2_C_ASC_Form )

    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources

    ! F => FA % Fluid_D ( )
    ! select type ( FF => F % Features )
    ! class is ( FluidFeatures_P_Form )
    ! call S % SetUseLimiter ( FF % Value ( :, FF % SHOCK ) )
    ! end select !-- FF

    end select !-- S


    !-- Template

    FinishTime = 1.0_KDR * TimeUnit
    if ( present ( FinishTimeOption ) ) &
      FinishTime = FinishTimeOption

    call FC % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, &
             FinishTimeOption = FinishTime, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

    call Show ( FC % Dimensionless, 'Dimensionless' )


    !-- Diagnostics
    
    call InitializeDiagnostics ( FC )

    
    !-- Cleanup

    end select !-- FA
    end select !-- GA
    end select !-- PS
    nullify ( F )

  end subroutine InitializeTemplate_FC


  subroutine OpenGridImageStreams ( I )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      I

    call I % OpenGridImageStreamsTemplate ( )

    if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return

    allocate ( I % GridImageStream_SA_EB )
    associate &
      ( GIS_I     => I % GridImageStream, &
        GIS_SA_EB => I % GridImageStream_SA_EB )
    call GIS_SA_EB % Initialize &
           ( trim ( I % Name ) // '_SA_EB', & 
             WorkingDirectoryOption = GIS_I % WorkingDirectory )
    end associate !-- GIS_I, etc.

  end subroutine OpenGridImageStreams


  subroutine OpenManifoldStreams ( I )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      I

    logical ( KDL ) :: &
      VerboseStream

    call I % OpenManifoldStreamsTemplate ( )

    if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return

    VerboseStream = .false.
    call PROGRAM_HEADER % GetParameter ( VerboseStream, 'VerboseStream' )

    associate ( GIS => I % GridImageStream_SA_EB )
    associate ( iS => 1 )  !-- iStream
      call I % PositionSpace_SA % OpenStream &
             ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
      call I % PositionSpace_EB % OpenStream &
             ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
    end associate !-- iS
    end associate !-- GIS

  end subroutine OpenManifoldStreams


  subroutine Write ( I )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      I

    call I % WriteTemplate ( )

    call ComputeSphericalAverage ( I )
    call ComputeEnclosedBaryons ( I )

    if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return

    associate &
      ( GIS => I % GridImageStream_SA_EB, &
        iS  => 1 )  !-- iStream
    call GIS % Open ( GIS % ACCESS_CREATE )

    call I % PositionSpace_SA % Write &
           ( iStream = iS, DirectoryOption = 'Chart_SA', &
             TimeOption = I % Time / I % TimeUnit, &
             CycleNumberOption = I % iCycle )
    call I % PositionSpace_EB % Write &
           ( iStream = iS, DirectoryOption = 'Chart_EB', &
             TimeOption = I % Time / I % TimeUnit, &
             CycleNumberOption = I % iCycle )

    call GIS % Close ( )
    end associate !-- GIS, etc.

  end subroutine Write


  subroutine FinalizeTemplate_FC ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    if ( allocated ( FC % Fluid_ASC_EB ) ) &
      deallocate ( FC % Fluid_ASC_EB )
    if ( allocated ( FC % Fluid_ASC_SA ) ) &
      deallocate ( FC % Fluid_ASC_SA )

    if ( allocated ( FC % PositionSpace_EB ) ) &
      deallocate ( FC % PositionSpace_EB )
    if ( allocated ( FC % PositionSpace_SA ) ) &
      deallocate ( FC % PositionSpace_SA )

    if ( allocated ( FC % GridImageStream_SA_EB ) ) &
      deallocate ( FC % GridImageStream_SA_EB )

    call FC % FinalizeTemplate_C_PS ( )

  end subroutine FinalizeTemplate_FC


  subroutine InitializeDiagnostics ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    type ( Real_1D_Form ), dimension ( 1 ) :: &
      Edge

    call Show ( 'Initializing FluidCentralCore diagnostics', &
                FC % IGNORABILITY )

    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( FA => FC % Current_ASC )
    class is ( Fluid_ASC_Form )


    !-- Spherical average position space

    allocate ( FC % PositionSpace_SA )
    associate ( PS_SA => FC % PositionSpace_SA )

    call PS_SA % Initialize ( 'SphericalAverage', nDimensionsOption = 1 )

    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )
      if ( FC % Dimensionless ) then
        call PS_SA % CreateChart &
               ( CoordinateSystemOption = 'SPHERICAL', &
                 nCellsOption           = C % nCells ( 1 : 1 ), &
                 nGhostLayersOption     = [ 0 ] )
      else
        call PS_SA % CreateChart &
               ( CoordinateSystemOption = 'SPHERICAL', &
                 CoordinateUnitOption   = FC % CoordinateUnit ( 1 : 1 ), &
                 nCellsOption           = C % nCells ( 1 : 1 ), &
                 nGhostLayersOption     = [ 0 ] )
      end if

      call Edge ( 1 ) % Initialize &
             ( C % Edge ( 1 ) % Value ( 1 : C % nCells ( 1 ) + 1 ) )

    end select !-- C

    call PS_SA % SetGeometry ( EdgeOption = Edge )


    !-- Enclosed baryons position space

    allocate ( FC % PositionSpace_EB )
    associate ( PS_EB => FC % PositionSpace_EB )

    call PS_EB % Initialize ( 'EnclosedBaryons', nDimensionsOption = 1 )

    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )
      if ( FC % Dimensionless ) then
        call PS_EB % CreateChart &
               ( CoordinateLabelOption = [ 'BaryonNumber' ], &
                 nCellsOption          = C % nCells ( 1 : 1 ), &
                 nGhostLayersOption    = [ 0 ] )
      else
        call PS_EB % CreateChart &
               ( CoordinateLabelOption = [ 'BaryonNumber' ], &
                 CoordinateUnitOption = [ UNIT % SOLAR_BARYON_NUMBER ], &
                 nCellsOption         = C % nCells ( 1 : 1 ), &
                 nGhostLayersOption   = [ 0 ] )
      end if
    end select !-- C

    call PS_EB % SetGeometry ( )


    !-- Spherical average fluid

    allocate ( FC % Fluid_ASC_SA )
    associate ( FA_SA => FC % Fluid_ASC_SA )

    if ( FC % Dimensionless ) then
      call FA_SA % Initialize ( PS_SA, FA % FluidType )
    else
      call FA_SA % Initialize &
             ( PS_SA, FA % FluidType, &
               AllocateTallyOption           =  .false., &
               AllocateSourcesOption         =  .false., &
               Velocity_U_UnitOption         =  FA % Velocity_U_Unit, &
               MomentumDensity_D_UnitOption  =  FA % MomentumDensity_D_Unit, &
               BaryonMassUnitOption          =  FA % BaryonMassUnit, &
               NumberDensityUnitOption       =  FA % NumberDensityUnit, &
               EnergyDensityUnitOption       =  FA % EnergyDensityUnit, &
               TemperatureUnitOption         =  FA % TemperatureUnit, &
               BaryonMassReferenceOption     =  FA % BaryonMassReference )
    end if


    !-- Enclosed baryons fluid

    allocate ( FC % Fluid_ASC_EB )
    associate ( FA_EB => FC % Fluid_ASC_EB )

    if ( FC % Dimensionless ) then
      call FA_EB % Initialize ( PS_EB, FA % FluidType )
    else
      call FA_EB % Initialize &
             ( PS_EB, FA % FluidType, &
               AllocateTallyOption           =  .false., &
               AllocateSourcesOption         =  .false., &
               Velocity_U_UnitOption         =  FA % Velocity_U_Unit, &
               MomentumDensity_D_UnitOption  =  FA % MomentumDensity_D_Unit, &
               BaryonMassUnitOption          =  FA % BaryonMassUnit, &
               NumberDensityUnitOption       =  FA % NumberDensityUnit, &
               EnergyDensityUnitOption       =  FA % EnergyDensityUnit, &
               BaryonMassReferenceOption     =  CONSTANT % ATOMIC_MASS_UNIT )
    end if


    !-- Cleanup

    end associate !-- FA_EB
    end associate !-- FA_SA
    end associate !-- PS_EB
    end associate !-- PS_SA
    end select !-- FA
    end select !-- PS

  end subroutine InitializeDiagnostics


  subroutine ComputeSphericalAverage ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    integer ( KDI ) :: &
      nAverage
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage
    type ( StorageForm ) :: &
      Integrand, &
      Average
    class ( GeometryFlatForm ), pointer :: &
      G_SA
    type ( SphericalAverageForm ) :: &
      SA
    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_SA

    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( FA => FC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )

    select type ( C_SA => FC % PositionSpace_SA % Chart )
    type is ( Chart_SLL_Form )

    G_SA => FC % PositionSpace_SA % Geometry ( )
    F_SA => FC % Fluid_ASC_SA % Fluid_D ( )

    nAverage = F % N_CONSERVED
    select type ( F )
    class is ( Fluid_P_HN_Form )
      nAverage = nAverage + 2  !-- Add pressure, temperature for initial guess
    end select !-- F

    allocate ( iaAverage ( nAverage ) )
    iaAverage ( 1 : F % N_CONSERVED ) = F % iaConserved
    select type ( F )
    class is ( Fluid_P_HN_Form )
      iaAverage ( F % N_CONSERVED + 1 )  =  F % PRESSURE
      iaAverage ( F % N_CONSERVED + 2 )  =  F % TEMPERATURE
    end select !-- F

    call Integrand % Initialize ( F, iaSelectedOption = iaAverage )
    call Average % Initialize ( F_SA, iaSelectedOption = iaAverage )

    call SA % Compute ( Average, C, C_SA, Integrand )

    if ( PS % Communicator % Rank == CONSOLE % DisplayRank ) &
      call F_SA % ComputeFromConserved ( G_SA )

    end select !-- C_SA
    end select !-- C
    end select !-- FA
    end select !-- PS
    nullify ( G_SA, F, F_SA )

  end subroutine ComputeSphericalAverage


  subroutine ComputeEnclosedBaryons ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    type ( Real_1D_Form ), dimension ( 1 ) :: &
      Edge
    class ( GeometryFlatForm ), pointer :: &
      G_SA
    class ( Fluid_D_Form ), pointer :: &
      F_SA, &
      F_EB

    integer ( KDI ) :: &
      iV

    if ( FC % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return

    G_SA => FC % PositionSpace_SA % Geometry ( )
    F_SA => FC % Fluid_ASC_SA % Fluid_D ( )
    F_EB => FC % Fluid_ASC_EB % Fluid_D ( )

    call Copy ( F_SA % Value, F_EB % Value )

    associate &
      ( dV => G_SA % Value ( :, G_SA % VOLUME ), &
         N => F_SA % Value ( :, F_SA % COMOVING_BARYON_DENSITY ) )

    call Edge ( 1 ) % Initialize ( size ( N ) + 1 )

    Edge ( 1 ) % Value ( 1 )  =  0.0_KDR
    do iV = 1, size ( N )
      Edge ( 1 ) % Value ( iV + 1 )  &
        =  Edge ( 1 ) % Value ( iV )  +  N ( iV ) * dV ( iV )
    end do

    select type ( C_EB => FC % PositionSpace_EB % Chart )
    type is ( Chart_SLL_Form )
      call C_EB % ResetGeometry ( Edge )
    end select !-- C_EB

    end associate !-- dV, etc.
    nullify ( G_SA, F_SA, F_EB )

  end subroutine ComputeEnclosedBaryons


  subroutine ApplySources &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( StorageForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call ApplyCurvilinear_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    call ApplyGravity_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

  end subroutine ApplySources


end module FluidCentral_Template
