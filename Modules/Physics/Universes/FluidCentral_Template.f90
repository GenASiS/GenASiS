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
        GridImageStream_SA_EM
      type ( Atlas_SC_Form ), allocatable :: &
        PositionSpace_SA, &  !-- SphericalAverage
        PositionSpace_EM     !-- EnclosedMass
      type ( Fluid_ASC_Form ), allocatable :: &
        Fluid_ASC_SA, &  !-- SphericalAverage
        Fluid_ASC_EM     !-- EnclosedMass
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
      ApplySources

      private :: &
        ComputeSolidAngleKernel, &
        ComputeMyIntegralKernel

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
      EnergyDensityUnit
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

    allocate ( I % GridImageStream_SA_EM )
    associate &
      ( GIS_I     => I % GridImageStream, &
        GIS_SA_EM => I % GridImageStream_SA_EM )
    call GIS_SA_EM % Initialize &
           ( trim ( I % Name ) // '_SA_EM', & 
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

    associate ( GIS => I % GridImageStream_SA_EM )
    associate ( iS => 1 )  !-- iStream
      call I % PositionSpace_SA % OpenStream &
             ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
      call I % PositionSpace_EM % OpenStream &
             ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
    end associate !-- iS
    end associate !-- GIS

  end subroutine OpenManifoldStreams


  subroutine Write ( I )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      I

    call I % WriteTemplate ( )

    call ComputeSphericalAverage ( I )

    if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return

    associate &
      ( GIS => I % GridImageStream_SA_EM, &
        iS  => 1 )  !-- iStream
    call GIS % Open ( GIS % ACCESS_CREATE )

    call I % PositionSpace_SA % Write &
           ( iStream = iS, DirectoryOption = 'Chart_SA', &
             TimeOption = I % Time / I % TimeUnit, &
             CycleNumberOption = I % iCycle )
    call I % PositionSpace_EM % Write &
           ( iStream = iS, DirectoryOption = 'Chart_EM', &
             TimeOption = I % Time / I % TimeUnit, &
             CycleNumberOption = I % iCycle )

    call GIS % Close ( )
    end associate !-- GIS, etc.

  end subroutine Write


  subroutine FinalizeTemplate_FC ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    if ( allocated ( FC % Fluid_ASC_EM ) ) &
      deallocate ( FC % Fluid_ASC_EM )
    if ( allocated ( FC % Fluid_ASC_SA ) ) &
      deallocate ( FC % Fluid_ASC_SA )

    if ( allocated ( FC % PositionSpace_EM ) ) &
      deallocate ( FC % PositionSpace_EM )
    if ( allocated ( FC % PositionSpace_SA ) ) &
      deallocate ( FC % PositionSpace_SA )

    if ( allocated ( FC % GridImageStream_SA_EM ) ) &
      deallocate ( FC % GridImageStream_SA_EM )

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

!    if ( PS % Communicator % Rank /= CONSOLE % DisplayRank ) &
!      return


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

      Edge ( 1 ) % Value  &
        =  C % Edge ( 1 ) % Value ( 1 : C % nCells ( 1 ) + 1 )

    end select !-- C

    call PS_SA % SetGeometry ( EdgeOption = Edge )


    !-- Enclosed mass position space

    allocate ( FC % PositionSpace_EM )
    associate ( PS_EM => FC % PositionSpace_EM )

    call PS_EM % Initialize ( 'EnclosedMass', nDimensionsOption = 1 )

    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )
      if ( FC % Dimensionless ) then
        call PS_EM % CreateChart &
               ( nCellsOption         = C % nCells ( 1 : 1 ), &
                 nGhostLayersOption   = C % nGhostLayers ( 1 : 1 ) )
      else
        call PS_EM % CreateChart &
               ( CoordinateUnitOption = [ UNIT % SOLAR_MASS ], &
                 nCellsOption         = C % nCells ( 1 : 1 ), &
                 nGhostLayersOption   = C % nGhostLayers ( 1 : 1 ) )
      end if
    end select !-- C

    call PS_EM % SetGeometry ( )


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
               BaryonMassReferenceOption     =  CONSTANT % ATOMIC_MASS_UNIT )
    end if


    !-- Enclosed mass fluid

    allocate ( FC % Fluid_ASC_EM )
    associate ( FA_EM => FC % Fluid_ASC_EM )

    if ( FC % Dimensionless ) then
      call FA_EM % Initialize ( PS_EM, FA % FluidType )
    else
      call FA_EM % Initialize &
             ( PS_EM, FA % FluidType, &
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

    end associate !-- FA_EM
    end associate !-- FA_SA
    end associate !-- PS_EM
    end associate !-- PS_SA
    end select !-- FA
    end select !-- PS

  end subroutine InitializeDiagnostics


  subroutine ComputeSphericalAverage ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    integer ( KDI ), dimension ( : ), allocatable :: &
      iRadius
    real ( KDR ), dimension ( : ), allocatable :: &
      SolidAngle
    real ( KDR ), dimension ( :, : ), allocatable :: &
      MyIntegral, &
      Integral
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_SA

    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( FA => FC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    select type ( C_SA => FC % PositionSpace_SA % Chart )
    type is ( Chart_SLL_Form )

    F_SA => FC % Fluid_ASC_SA % Fluid_D ( )

    allocate ( iRadius ( G % nValues ) )
    allocate ( SolidAngle ( G % nValues ) )
    allocate ( MyIntegral ( C_SA % nCells ( 1 ), F_SA % N_CONSERVED ) )
    allocate ( Integral ( C_SA % nCells ( 1 ), F_SA % N_CONSERVED ) )

    call ComputeSolidAngleKernel &
           ( SolidAngle, iRadius, &
             G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), &
             G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
             G % Value ( :, G % CENTER_U ( 1 ) ), &
             G % Value ( :, G % CENTER_U ( 2 ) ), &
             C_SA % Edge ( 1 ) % Value, PS % nDimensions )

    call Clear ( MyIntegral )

    end select !-- C_SA
    end select !-- FA
    end select !-- PS
    nullify ( G, F, F_SA )

  end subroutine ComputeSphericalAverage


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


  subroutine ComputeSolidAngleKernel &
               ( dOmega, iR, W_L_2, W_L_3, W_R_2, W_R_3, R_C, Th_C, &
                 EdgeValue, nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dOmega
    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      iR
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      W_L_2, W_L_3, &
      W_R_2, W_R_3, &
      R_C, &
      Th_C, &
      EdgeValue
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      TwoPi, FourPi, &
      Th_I, Th_O, &
      dPh
    
    TwoPi   =  2.0_KDR  *  CONSTANT % PI
    FourPi  =  4.0_KDR  *  CONSTANT % PI

    !$OMP parallel do private ( iV, Th_I, Th_O, dPh )
    do iV = 1, size ( dOmega )
      
      Th_I  =  Th_C ( iV )  -  W_L_2 ( iV )
      Th_O  =  Th_C ( iV )  +  W_R_2 ( iV )

      dPh  =  W_L_3 ( iV )  +  W_R_3 ( iV ) 

      select case ( nDimensions )
      case ( 1 )
        dOmega ( iV )  =  FourPi
      case ( 2 )
        dOmega ( iV )  =  TwoPi  *  ( cos ( Th_I )  -  cos ( Th_O ) )
      case ( 3 ) 
        dOmega ( iV )  =  dPh  *  ( cos ( Th_I )  -  cos ( Th_O ) )
      end select !-- nDimensions

      call Search ( EdgeValue, R_C ( iV ), iR ( iV ) )

    end do
    !$OMP end parallel do

  end subroutine ComputeSolidAngleKernel


  subroutine ComputeMyIntegralKernel ( MyI, IsProperCell, D, dOmega, iR )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      MyI
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      D, &
      dOmega
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iR

  end subroutine ComputeMyIntegralKernel


end module FluidCentral_Template
