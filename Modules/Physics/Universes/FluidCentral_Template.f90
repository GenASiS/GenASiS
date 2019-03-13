module FluidCentral_Template

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ComputeGravity_Command
  use ApplyGravity_F__Command
  use Universe_Template

  implicit none
  private

  type, public, extends ( UniverseTemplate ), abstract :: &
    FluidCentralTemplate
!       integer ( KDI ), dimension ( : ), allocatable :: &
!         nCoarsen_2, &
!         nCoarsen_3
      real ( KDR ) :: &
        RadiusPolarMomentum = 0.0_KDR
      logical ( KDL ) :: &
        Dimensionless, &
        UseCoarsening
!       type ( StorageForm ), dimension ( : ), allocatable :: &
!         CoarsenPillar_2, &
!         CoarsenPillar_3
!       type ( CollectiveOperation_R_Form ), allocatable :: &
!         CO_CoarsenForward_2, CO_CoarsenBackward_2, &
!         CO_CoarsenForward_3, CO_CoarsenBackward_3
!       type ( GridImageStreamForm ), allocatable :: &
!         GridImageStream_SA_EB
!       type ( Atlas_SC_Form ), allocatable :: &
!         PositionSpace_SA, &  !-- SphericalAverage
!         PositionSpace_EB     !-- EnclosedBaryons
!       type ( Fluid_ASC_Form ), allocatable :: &
!         Fluid_ASC_SA, &  !-- SphericalAverage
!         Fluid_ASC_EB     !-- EnclosedBaryons
  contains
    procedure, public, pass :: &
      InitializeTemplate_FC
    procedure, public, pass :: &
      FinalizeTemplate_FC
    procedure, private, pass :: &
      AllocateIntegrator_FC
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_FC
    procedure, public, pass :: &
      InitializePositionSpace
    procedure ( IA ), private, pass, deferred :: &
      InitializeAtlas
    procedure ( IG ), private, pass, deferred :: &
      InitializeGeometry
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      InitializeStep
!     procedure, public, pass :: &
!       SetCoarseningTemplate
!     procedure ( SC ), private, pass, deferred :: &
!       SetCoarsening
!     procedure, public, pass :: &
!       CoarsenSingularityTemplate
!     procedure ( CS ), public, nopass, deferred :: &
!       CoarsenSingularities
!     procedure, public, pass :: &  !-- 2
!       OpenGridImageStreams
!     procedure, public, pass :: &  !-- 2
!       OpenManifoldStreams
!     procedure, public, pass :: &  !-- 3
!       Write
  end type FluidCentralTemplate

  abstract interface

    subroutine IA ( FC, RadiusMaxOption, RadiusCoreOption, &
                    RadiusExcisionOption, RadialRatioOption, &
                    nCellsPolarOption )
      use Basics
      import FluidCentralTemplate
      class ( FluidCentralTemplate ), intent ( inout ) :: &
        FC
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusExcisionOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsPolarOption
    end subroutine IA

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

!     subroutine SC ( FC )
!       import FluidCentralTemplate
!       class ( FluidCentralTemplate ), intent ( inout ) :: &
!         FC
!     end subroutine SC

!     subroutine CS ( S )
!       use Basics
!       class ( StorageForm ), intent ( inout ) :: &
!         S
!     end subroutine CS

  end interface

    private :: &
!       InitializeDiagnostics, &
!       ComputeSphericalAverage, &
!       ComputeEnclosedBaryons, &
!       ComposePillars, &
!       CoarsenPillars, &
!       DecomposePillars, &
      ApplySources

contains


  subroutine InitializeTemplate_FC &
               ( FC, Name, FluidType, GeometryType, DimensionlessOption, &
                 FinishTimeOption, CourantFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, & 
                 RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
                 RadialRatioOption, CentralMassOption, &
                 nWriteOption, nCellsPolarOption )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
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
      nCellsPolarOption

    real ( KDR ) :: &
      FinishTime

    if ( FC % Type == '' ) &
      FC % Type = 'a FluidCentral'

    FC % Dimensionless = .false.
    if ( present ( DimensionlessOption ) ) &
      FC % Dimensionless = DimensionlessOption

    call FC % InitializeTemplate ( Name )

    call FC % AllocateIntegrator &
           ( )
    call FC % InitializePositionSpace &
           ( GeometryType, RadiusMaxOption, RadiusCoreOption, &
             RadiusMinOption, RadialRatioOption, CentralMassOption, &
             nCellsPolarOption )
    call FC % InitializeFluid &
           ( FluidType, LimiterParameterOption, ShockThresholdOption )
    call FC % InitializeStep &
           ( Name ) 

    FinishTime  =  1.0_KDR  *  FC % Units % Time
    if ( present ( FinishTimeOption ) ) &
      FinishTime = FinishTimeOption

    select type ( I => FC % Integrator )
    class is ( Integrator_C_PS_Form )
      call I % Initialize &
             ( Name, TimeUnitOption = FC % Units % Time, &
               FinishTimeOption = FinishTime, &
               CourantFactorOption = CourantFactorOption, &
               nWriteOption = nWriteOption )
    end select !-- I

    call Show ( FC % Dimensionless, 'Dimensionless', FC % IGNORABILITY )
    call Show ( FC % UseCoarsening, 'UseCoarsening', FC % IGNORABILITY )
    call Show ( FC % RadiusPolarMomentum, FC % Units % Coordinate_PS ( 1 ), &
                'RadiusPolarMomentum', FC % IGNORABILITY )

!     if ( FC % UseCoarsening ) &
!       call FC % SetCoarsening ( )


!     !-- Diagnostics
    
!     call InitializeDiagnostics ( FC )

    
!     !-- Cleanup

!     end select !-- FA
!     end select !-- GA
!     end select !-- PS
!     nullify ( F )

  end subroutine InitializeTemplate_FC


  subroutine FinalizeTemplate_FC ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

!     if ( allocated ( FC % Fluid_ASC_EB ) ) &
!       deallocate ( FC % Fluid_ASC_EB )
!     if ( allocated ( FC % Fluid_ASC_SA ) ) &
!       deallocate ( FC % Fluid_ASC_SA )

!     if ( allocated ( FC % PositionSpace_EB ) ) &
!       deallocate ( FC % PositionSpace_EB )
!     if ( allocated ( FC % PositionSpace_SA ) ) &
!       deallocate ( FC % PositionSpace_SA )

!     if ( allocated ( FC % GridImageStream_SA_EB ) ) &
!       deallocate ( FC % GridImageStream_SA_EB )

!     if ( allocated ( FC % CO_CoarsenBackward_3 ) ) &
!       deallocate ( FC % CO_CoarsenBackward_3 )
!     if ( allocated ( FC % CO_CoarsenForward_3 ) ) &
!       deallocate ( FC % CO_CoarsenForward_3 )
!     if ( allocated ( FC % CO_CoarsenBackward_2 ) ) &
!       deallocate ( FC % CO_CoarsenBackward_2 )
!     if ( allocated ( FC % CO_CoarsenForward_2 ) ) &
!       deallocate ( FC % CO_CoarsenForward_2 )

!     if ( allocated ( FC % CoarsenPillar_3 ) ) &
!       deallocate ( FC % CoarsenPillar_3 )
!     if ( allocated ( FC % CoarsenPillar_2 ) ) &
!       deallocate ( FC % CoarsenPillar_2 )

!     if ( allocated ( FC % nCoarsen_3 ) ) &
!       deallocate ( FC % nCoarsen_3 )
!     if ( allocated ( FC % nCoarsen_2 ) ) &
!       deallocate ( FC % nCoarsen_2 )

    call FC % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_FC


  subroutine AllocateIntegrator_FC ( FC )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC

    allocate ( Integrator_C_PS_Form :: FC % Integrator )

  end subroutine AllocateIntegrator_FC


  subroutine InitializePositionSpace &
               ( FC, GeometryType, RadiusMaxOption, RadiusCoreOption, &
                 RadiusMinOption, RadialRatioOption, CentralMassOption, &
                 nCellsPolarOption )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in )  :: &
      GeometryType
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusMinOption, &
      RadialRatioOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    if ( .not. FC % Dimensionless ) then
      FC % Units % Coordinate_PS  &
        =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]
      FC % Units % Time &
        =  UNIT % SECOND
    end if

    call FC % InitializeAtlas &
           ( RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
             RadialRatioOption, nCellsPolarOption )

    select type ( PS => FC % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    call FC % InitializeGeometry ( GA, PS, GeometryType, CentralMassOption )

    call PS % SetGeometry ( GA )

    FC % UseCoarsening = .true.
    call PROGRAM_HEADER % GetParameter ( FC % UseCoarsening, 'UseCoarsening' )
    if ( FC % UseCoarsening ) &
      call PS % SetCoarsening ( )

    end select !-- GA
    end select !-- PS

  end subroutine InitializePositionSpace


  subroutine InitializeFluid &
               ( FC, FluidType, LimiterParameterOption, ShockThresholdOption )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in )  :: &
      FluidType
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption, &
      ShockThresholdOption

    real ( KDR ) :: &
      BaryonMassReference

    select type ( I => FC % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Fluid_ASC_Form :: I % Current_ASC )
    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    if ( .not. FC % Dimensionless ) then

      FC % Units % BaryonMass     =  UNIT % ATOMIC_MASS_UNIT
      FC % Units % NumberDensity  =  UNIT % NUMBER_DENSITY_NUCLEAR
      FC % Units % EnergyDensity  =  UNIT % ENERGY_DENSITY_NUCLEAR
      FC % Units % Temperature    =  UNIT % MEGA_ELECTRON_VOLT

      FC % Units % Velocity_U  &
        =  FC % Units % Coordinate_PS  /  FC % Units % Time

      FC % Units % MomentumDensity_D  &
        =  FC % Units % BaryonMass  *  FC % Units % NumberDensity  &
           *  FC % Units % Velocity_U
      FC % Units % MomentumDensity_D ( 2 )  &
        =  FC % Units % MomentumDensity_D ( 2 )  &
           *  FC % Units % Coordinate_PS ( 1 ) ** 2
      FC % Units % MomentumDensity_D ( 3 )  &
        =  FC % Units % MomentumDensity_D ( 3 )  &
           *  FC % Units % Coordinate_PS ( 1 ) ** 2

      FC % Units % MomentumDensity_U  &
        =  FC % Units % BaryonMass  *  FC % Units % NumberDensity  &
           *  FC % Units % Velocity_U
      FC % Units % MomentumDensity_U ( 2 )  &
        =  FC % Units % MomentumDensity_U ( 2 )  &
           /  FC % Units % Coordinate_PS ( 1 ) ** 2
      FC % Units % MomentumDensity_U ( 3 )  &
        =  FC % Units % MomentumDensity_U ( 3 )  &
           /  FC % Units % Coordinate_PS ( 1 ) ** 2

      FC % Units % Number           =  UNIT % SOLAR_BARYON_NUMBER
      FC % Units % Energy           =  UNIT % ENERGY_SOLAR_MASS
      FC % Units % Momentum         =  UNIT % MOMENTUM_SOLAR_MASS
      FC % Units % AngularMomentum  =  UNIT % SOLAR_KERR_PARAMETER

      BaryonMassReference  =  CONSTANT % ATOMIC_MASS_UNIT

    else
      BaryonMassReference  =  1.0_KDR
    end if

    call FA % Initialize &
           ( PS, FluidType, FC % Units, &
             BaryonMassReferenceOption = BaryonMassReference, &
             LimiterParameterOption = LimiterParameterOption, &
             ShockThresholdOption = ShockThresholdOption )

    end select !-- FA
    end select !-- PS
    end select !-- I

  end subroutine InitializeFluid


  subroutine InitializeStep ( FC, Name )

    class ( FluidCentralTemplate ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in )  :: &
      Name

    select type ( I => FC % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    allocate ( Step_RK2_C_ASC_Form :: I % Step )
    select type ( S => I % Step )
    class is ( Step_RK2_C_ASC_Form )

    call S % Initialize ( I, FA, Name )
    S % ComputeConstraints % Pointer => ComputeGravity
    S % ApplySources % Pointer => ApplySources

    ! F => FA % Fluid_D ( )
    ! select type ( FF => F % Features )
    ! class is ( FluidFeatures_P_Form )
    ! call S % SetUseLimiter ( FF % Value ( :, FF % SHOCK ) )
    ! end select !-- FF

    end select !-- S
    end select !-- FA
    end select !-- I

  end subroutine InitializeStep


!   subroutine SetCoarseningTemplate ( FC, iAngular )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC
!     integer ( KDI ), intent ( in ) :: &
!       iAngular

!     integer ( KDI ) :: &
!       iP, &  !-- iPillar
!       nValuesFactor_F_2, nValuesFactor_B_2, &
!       nValuesFactor_F_3, nValuesFactor_B_3
!     class ( Fluid_D_Form ), pointer :: &
!       F

!     select type ( FA => FC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_D ( )

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_C_Template )

!     select case ( iAngular )
!     case ( 2 )

!       if ( .not. C % Communicator_2 % Initialized ) &
!         return

!       allocate ( FC % CO_CoarsenForward_2 )
!       allocate ( FC % CO_CoarsenBackward_2 )
!       associate &
!         ( CO_F => FC % CO_CoarsenForward_2, &
!           CO_B => FC % CO_CoarsenBackward_2 )

!       nValuesFactor_F_2  =  ( 1 + F % N_CONSERVED )  *  C % nCellsBrick ( 2 ) &
!                             +  1
!       nValuesFactor_B_2  =  F % N_CONSERVED  *  C % nCellsBrick ( 2 )

!       call CO_F % Initialize &
!              ( C % Communicator_2, &
!                nIncoming  =  C % nSegmentsFrom_2  *  nValuesFactor_F_2, &
!                nOutgoing  =  C % nSegmentsTo_2    *  nValuesFactor_F_2 )
!       call CO_B % Initialize &
!              ( C % Communicator_2, &
!                nIncoming  =  C % nSegmentsTo_2    *  nValuesFactor_B_2, &
!                nOutgoing  =  C % nSegmentsFrom_2  *  nValuesFactor_B_2 )

!       end associate !-- CO_F, etc.

!       allocate ( FC % nCoarsen_2 ( C % nPillars_2 ) )
!       allocate ( FC % CoarsenPillar_2 ( C % nPillars_2 ) )
!       do iP = 1, C % nPillars_2
!         call FC % CoarsenPillar_2 ( iP ) % Initialize &
!                ( [ C % nCells ( 2 ), 1 + F % N_CONSERVED ] )
!       end do  !-- iP

!     case ( 3 )

!       if ( .not. C % Communicator_3 % Initialized ) &
!         return

!       allocate ( FC % CO_CoarsenForward_3 )
!       allocate ( FC % CO_CoarsenBackward_3 )
!       associate &
!         ( CO_F => FC % CO_CoarsenForward_3, &
!           CO_B => FC % CO_CoarsenBackward_3 )

!       nValuesFactor_F_3  =  ( 1 + F % N_CONSERVED )  *  C % nCellsBrick ( 3 ) &
!                             +  1
!       nValuesFactor_B_3  =  F % N_CONSERVED  *  C % nCellsBrick ( 3 )
  
!       call CO_F % Initialize &
!              ( C % Communicator_3, &
!                nIncoming  =  C % nSegmentsFrom_3  *  nValuesFactor_F_3, &
!                nOutgoing  =  C % nSegmentsTo_3    *  nValuesFactor_F_3 )
!       call CO_B % Initialize &
!              ( C % Communicator_3, &
!                nIncoming  =  C % nSegmentsTo_3    *  nValuesFactor_B_3, &
!                nOutgoing  =  C % nSegmentsFrom_3  *  nValuesFactor_B_3 )

!       end associate !-- CO_F, etc.

!       allocate ( FC % nCoarsen_3 ( C % nPillars_3 ) )
!       allocate ( FC % CoarsenPillar_3 ( C % nPillars_3 ) )
!       do iP = 1, C % nPillars_3
!         call FC % CoarsenPillar_3 ( iP ) % Initialize &
!                ( [ C % nCells ( 3 ), 1 + F % N_CONSERVED ] )
!       end do  !-- iP

!     end select !-- iAngular
!     end select !-- C
!     end select !-- PS
!     end select !-- FA
!     nullify ( F )

!   end subroutine SetCoarseningTemplate


!   subroutine CoarsenSingularityTemplate ( FC, S, iAngular )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC
!     class ( StorageForm ), intent ( inout ) :: &
!       S
!     integer ( KDI ), intent ( in ) :: &
!       iAngular

!     if ( .not. FC % UseCoarsening ) &
!       return

!     call Show ( 'Coarsening Singularity', FC % IGNORABILITY + 2 )
!     call Show ( iAngular, 'iAngular', FC % IGNORABILITY + 2 )

!     call ComposePillars ( FC, S, iAngular )
!     call CoarsenPillars ( FC, iAngular )
!     call DecomposePillars ( FC, S, iAngular )

!   end subroutine CoarsenSingularityTemplate


!   subroutine OpenGridImageStreams ( I )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       I

!     call I % OpenGridImageStreamsTemplate ( )

!     if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
!       return

!     allocate ( I % GridImageStream_SA_EB )
!     associate &
!       ( GIS_I     => I % GridImageStream, &
!         GIS_SA_EB => I % GridImageStream_SA_EB )
!     call GIS_SA_EB % Initialize &
!            ( trim ( I % Name ) // '_SA_EB', & 
!              WorkingDirectoryOption = GIS_I % WorkingDirectory )
!     end associate !-- GIS_I, etc.

!   end subroutine OpenGridImageStreams


!   subroutine OpenManifoldStreams ( I, VerboseStreamOption )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       I
!     logical ( KDL ), intent ( in ), optional :: &
!       VerboseStreamOption

!     logical ( KDL ) :: &
!       VerboseStream

!     call I % OpenManifoldStreamsTemplate ( VerboseStreamOption )

!     if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
!       return

!     VerboseStream = .false.
!     if ( present ( VerboseStreamOption ) ) &
!       VerboseStream = VerboseStreamOption
!     call PROGRAM_HEADER % GetParameter ( VerboseStream, 'VerboseStream' )

!     associate ( GIS => I % GridImageStream_SA_EB )
!     associate ( iS => 1 )  !-- iStream
!       call I % PositionSpace_SA % OpenStream &
!              ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
!       call I % PositionSpace_EB % OpenStream &
!              ( GIS, 'Time', iStream = iS, VerboseOption = VerboseStream )
!     end associate !-- iS
!     end associate !-- GIS

!   end subroutine OpenManifoldStreams


!   subroutine Write ( I )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       I

!     call I % WriteTemplate ( )

!     call ComputeSphericalAverage ( I )
!     call ComputeEnclosedBaryons ( I )

!     if ( I % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
!       return

!     associate &
!       ( GIS => I % GridImageStream_SA_EB, &
!         iS  => 1 )  !-- iStream
!     call GIS % Open ( GIS % ACCESS_CREATE )

!     call I % PositionSpace_SA % Write &
!            ( iStream = iS, DirectoryOption = 'Chart_SA', &
!              TimeOption = I % Time / I % TimeUnit, &
!              CycleNumberOption = I % iCycle )
!     call I % PositionSpace_EB % Write &
!            ( iStream = iS, DirectoryOption = 'Chart_EB', &
!              TimeOption = I % Time / I % TimeUnit, &
!              CycleNumberOption = I % iCycle )

!     call GIS % Close ( )
!     end associate !-- GIS, etc.

!   end subroutine Write


!   subroutine InitializeDiagnostics ( FC )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC

!     type ( Real_1D_Form ), dimension ( 1 ) :: &
!       Edge

!     call Show ( 'Initializing FluidCentralCore diagnostics', &
!                 FC % IGNORABILITY )

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( FA => FC % Current_ASC )
!     class is ( Fluid_ASC_Form )


!     !-- Spherical average position space

!     allocate ( FC % PositionSpace_SA )
!     associate ( PS_SA => FC % PositionSpace_SA )

!     call PS_SA % Initialize ( 'SphericalAverage', nDimensionsOption = 1 )

!     select type ( C => PS % Chart )
!     class is ( Chart_SL_Template )
!       if ( FC % Dimensionless ) then
!         call PS_SA % CreateChart &
!                ( CoordinateSystemOption = 'SPHERICAL', &
!                  nCellsOption           = C % nCells ( 1 : 1 ), &
!                  nGhostLayersOption     = [ 0 ] )
!       else
!         call PS_SA % CreateChart &
!                ( CoordinateSystemOption = 'SPHERICAL', &
!                  CoordinateUnitOption   = FC % CoordinateUnit ( 1 : 1 ), &
!                  nCellsOption           = C % nCells ( 1 : 1 ), &
!                  nGhostLayersOption     = [ 0 ] )
!       end if

!       call Edge ( 1 ) % Initialize &
!              ( C % Edge ( 1 ) % Value ( 1 : C % nCells ( 1 ) + 1 ) )

!     end select !-- C

!     call PS_SA % SetGeometry ( EdgeOption = Edge )


!     !-- Enclosed baryons position space

!     allocate ( FC % PositionSpace_EB )
!     associate ( PS_EB => FC % PositionSpace_EB )

!     call PS_EB % Initialize ( 'EnclosedBaryons', nDimensionsOption = 1 )

!     select type ( C => PS % Chart )
!     class is ( Chart_SL_Template )
!       if ( FC % Dimensionless ) then
!         call PS_EB % CreateChart &
!                ( CoordinateLabelOption = [ 'BaryonNumber' ], &
!                  nCellsOption          = C % nCells ( 1 : 1 ), &
!                  nGhostLayersOption    = [ 0 ] )
!       else
!         call PS_EB % CreateChart &
!                ( CoordinateLabelOption = [ 'BaryonNumber' ], &
!                  CoordinateUnitOption = [ UNIT % SOLAR_BARYON_NUMBER ], &
!                  nCellsOption         = C % nCells ( 1 : 1 ), &
!                  nGhostLayersOption   = [ 0 ] )
!       end if
!     end select !-- C

!     call PS_EB % SetGeometry ( )


!     !-- Spherical average fluid

!     allocate ( FC % Fluid_ASC_SA )
!     associate ( FA_SA => FC % Fluid_ASC_SA )

!     if ( FC % Dimensionless ) then
!       call FA_SA % Initialize &
!              ( PS_SA, FA % FluidType, &
!                AllocateTallyOption           =  .false., &
!                AllocateSourcesOption         =  .false. )
!     else
!       call FA_SA % Initialize &
!              ( PS_SA, FA % FluidType, &
!                AllocateTallyOption           =  .false., &
!                AllocateSourcesOption         =  .false., &
!                Velocity_U_UnitOption         =  FA % Velocity_U_Unit, &
!                MomentumDensity_D_UnitOption  =  FA % MomentumDensity_D_Unit, &
!                BaryonMassUnitOption          =  FA % BaryonMassUnit, &
!                NumberDensityUnitOption       =  FA % NumberDensityUnit, &
!                EnergyDensityUnitOption       =  FA % EnergyDensityUnit, &
!                TemperatureUnitOption         =  FA % TemperatureUnit, &
!                BaryonMassReferenceOption     =  FA % BaryonMassReference )
!     end if


!     !-- Enclosed baryons fluid

!     allocate ( FC % Fluid_ASC_EB )
!     associate ( FA_EB => FC % Fluid_ASC_EB )

!     if ( FC % Dimensionless ) then
!       call FA_EB % Initialize &
!              ( PS_EB, FA % FluidType, &
!                AllocateTallyOption           =  .false., &
!                AllocateSourcesOption         =  .false. )
!     else
!       call FA_EB % Initialize &
!              ( PS_EB, FA % FluidType, &
!                AllocateTallyOption           =  .false., &
!                AllocateSourcesOption         =  .false., &
!                Velocity_U_UnitOption         =  FA % Velocity_U_Unit, &
!                MomentumDensity_D_UnitOption  =  FA % MomentumDensity_D_Unit, &
!                BaryonMassUnitOption          =  FA % BaryonMassUnit, &
!                NumberDensityUnitOption       =  FA % NumberDensityUnit, &
!                EnergyDensityUnitOption       =  FA % EnergyDensityUnit, &
!                BaryonMassReferenceOption     =  CONSTANT % ATOMIC_MASS_UNIT )
!     end if


!     !-- Cleanup

!     end associate !-- FA_EB
!     end associate !-- FA_SA
!     end associate !-- PS_EB
!     end associate !-- PS_SA
!     end select !-- FA
!     end select !-- PS

!   end subroutine InitializeDiagnostics


!   subroutine ComputeSphericalAverage ( FC )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC

!     integer ( KDI ) :: &
!       nAverage
!     integer ( KDI ), dimension ( : ), allocatable :: &
!       iaAverage
!     type ( StorageForm ) :: &
!       Integrand, &
!       Average
!     class ( GeometryFlatForm ), pointer :: &
!       G_SA
!     type ( SphericalAverageForm ) :: &
!       SA
!     class ( Fluid_D_Form ), pointer :: &
!       F, &
!       F_SA

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( FA => FC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_D ( )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_Form )

!     select type ( C_SA => FC % PositionSpace_SA % Chart )
!     type is ( Chart_SLL_Form )

!     G_SA => FC % PositionSpace_SA % Geometry ( )
!     F_SA => FC % Fluid_ASC_SA % Fluid_D ( )

!     nAverage = F % N_CONSERVED
!     select type ( F )
!     class is ( Fluid_P_HN_Form )
!       nAverage = nAverage + 2  !-- Add pressure, temperature for initial guess
!     end select !-- F

!     allocate ( iaAverage ( nAverage ) )
!     iaAverage ( 1 : F % N_CONSERVED ) = F % iaConserved
!     select type ( F )
!     class is ( Fluid_P_HN_Form )
!       iaAverage ( F % N_CONSERVED + 1 )  =  F % PRESSURE
!       iaAverage ( F % N_CONSERVED + 2 )  =  F % TEMPERATURE
!     end select !-- F

!     call Integrand % Initialize ( F, iaSelectedOption = iaAverage )
!     call Average % Initialize ( F_SA, iaSelectedOption = iaAverage )

!     call SA % Compute ( Average, C, C_SA, Integrand )

!     if ( PS % Communicator % Rank == CONSOLE % DisplayRank ) &
!       call F_SA % ComputeFromConserved ( G_SA )

!     end select !-- C_SA
!     end select !-- C
!     end select !-- FA
!     end select !-- PS
!     nullify ( G_SA, F, F_SA )

!   end subroutine ComputeSphericalAverage


!   subroutine ComputeEnclosedBaryons ( FC )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC

!     type ( Real_1D_Form ), dimension ( 1 ) :: &
!       Edge
!     class ( GeometryFlatForm ), pointer :: &
!       G_SA
!     class ( Fluid_D_Form ), pointer :: &
!       F_SA, &
!       F_EB

!     integer ( KDI ) :: &
!       iV

!     if ( FC % PositionSpace % Communicator % Rank /= CONSOLE % DisplayRank ) &
!       return

!     G_SA => FC % PositionSpace_SA % Geometry ( )
!     F_SA => FC % Fluid_ASC_SA % Fluid_D ( )
!     F_EB => FC % Fluid_ASC_EB % Fluid_D ( )

!     call Copy ( F_SA % Value, F_EB % Value )

!     associate &
!       ( dV => G_SA % Value ( :, G_SA % VOLUME ), &
!          N => F_SA % Value ( :, F_SA % COMOVING_BARYON_DENSITY ) )

!     call Edge ( 1 ) % Initialize ( size ( N ) + 1 )

!     Edge ( 1 ) % Value ( 1 )  =  0.0_KDR
!     do iV = 1, size ( N )
!       Edge ( 1 ) % Value ( iV + 1 )  &
!         =  Edge ( 1 ) % Value ( iV )  +  N ( iV ) * dV ( iV )
!     end do

!     select type ( C_EB => FC % PositionSpace_EB % Chart )
!     type is ( Chart_SLL_Form )
!       call C_EB % ResetGeometry ( Edge )
!     end select !-- C_EB

!     end associate !-- dV, etc.
!     nullify ( G_SA, F_SA, F_EB )

!   end subroutine ComputeEnclosedBaryons


!   subroutine ComposePillars ( FC, S, iAngular )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC
!     class ( StorageForm ), intent ( inout ) :: &
!       S
!     integer ( KDI ), intent ( in ) :: &
!       iAngular

!     integer ( KDI ) :: &
!       oO, oI, &      !-- oOutgoing, oIncoming
!       oP, &          !-- oPillar
!       oC, &          !-- oCell
!       iC, jC, kC, &  !-- iCell, etc.
!       iS, &          !-- iSelected
!       iG, &          !-- iGroup (of pillars/processes)
!       iB, &          !-- iBrick
!       iR, &          !-- iRank
!       iPS, &         !-- iPillarSegment
!       iP             !-- iPillar
!     real ( KDR ), dimension ( :, :, : ), pointer :: &
!       Crsn_2, Crsn_3, &
!       Vol, &
!       SV
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_Form )
!     G => PS % Geometry ( )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_C_Template )

!     associate ( nCB => C % nCellsBrick )

!     select case ( iAngular )
!     case ( 2 )

!       if ( .not. C % Communicator_2 % Initialized ) &
!         return

!       associate &
!         ( CO => FC % CO_CoarsenForward_2 )
!       associate &
!         ( Outgoing => CO % Outgoing % Value, &
!           Incoming => CO % Incoming % Value )

!       call C % SetVariablePointer &
!              ( G % Value ( :, G % COARSENING ( 2 ) ), Crsn_2 )
!       call C % SetVariablePointer &
!              ( G % Value ( :, G % VOLUME ), Vol )

!       oO = 0
!       do kC = 1, nCB ( 3 )
!         do iC = 1, nCB ( 1 )
!           if ( Crsn_2 ( iC, 1, kC ) < 2.0_KDR ) &
!             cycle
!           Outgoing ( oO + 1 ) = Crsn_2 ( iC, 1, kC )
!           oO = oO + 1
!           Outgoing ( oO + 1 : oO + nCB ( 2 ) ) &
!             =  Vol ( iC, 1 : nCB ( 2 ), kC )
!           oO = oO + nCB ( 2 )
!           do iS = 1, S % nVariables
!             call C % SetVariablePointer &
!                    ( S % Value ( :, S % iaSelected ( iS ) ), SV )
!             Outgoing ( oO + 1 : oO + nCB ( 2 ) ) &
!               =  SV ( iC, 1 : nCB ( 2 ), kC )
!             oO = oO + nCB ( 2 )
!           end do !-- iS
!         end do !-- iC
!       end do !-- kC

!       call CO % AllToAll_V ( )

!       oP = 0
!       oI = 0
!       do iG = 1, C % Communicator_2 % Size  /  C % nBricks ( 2 )
!         oC = 0
!         do iB = 1, C % nBricks ( 2 )
!           iR = ( iG - 1 ) * C % nBricks ( 2 )  +  iB  -  1
!           do iPS = 1, C % nSegmentsFrom_2 ( iR )
!             iP = oP + iPS
!             FC % nCoarsen_2 ( iP ) = int ( Incoming ( oI + 1 ) + 0.5_KDR )
!             oI = oI + 1
!             associate ( CP => FC % CoarsenPillar_2 ( iP ) )
!             do iS = 1, CP % nVariables
!               CP % Value ( oC + 1 : oC + nCB ( 2 ), iS )  &
!                 =  Incoming ( oI + 1 : oI + nCB ( 2 ) )
!               oI = oI + nCB ( 2 )
!             end do !-- iS
!             end associate !-- CP
!           end do !-- iPS
!           oC = oC + nCB ( 2 )
!         end do !-- iB
!         oP = oP + C % nSegmentsFrom_2 ( iR )
!       end do !-- iG

!       end associate !-- Outgoing, etc.
!       end associate !-- CO

!     case ( 3 )

!       if ( .not. C % Communicator_3 % Initialized ) &
!         return

!       associate &
!         ( CO => FC % CO_CoarsenForward_3 )
!       associate &
!         ( Outgoing => CO % Outgoing % Value, &
!           Incoming => CO % Incoming % Value )

!       call C % SetVariablePointer &
!              ( G % Value ( :, G % COARSENING ( 3 ) ), Crsn_3 )
!       call C % SetVariablePointer &
!              ( G % Value ( :, G % VOLUME ), Vol )

!       oO = 0
!       do jC = 1, nCB ( 2 )
!         do iC = 1, nCB ( 1 )
!           if ( Crsn_3 ( iC, jC, 1 ) < 2.0_KDR ) &
!             cycle
!           Outgoing ( oO + 1 ) = Crsn_3 ( iC, jC, 1 )
!           oO = oO + 1
!           Outgoing ( oO + 1 : oO + nCB ( 3 ) ) &
!             =  Vol ( iC, jC, 1 : nCB ( 3 ) )
!           oO = oO + nCB ( 3 )
!           do iS = 1, S % nVariables
!             call C % SetVariablePointer &
!                    ( S % Value ( :, S % iaSelected ( iS ) ), SV )
!             Outgoing ( oO + 1 : oO + nCB ( 3 ) ) &
!               =  SV ( iC, jC, 1 : nCB ( 3 ) )
!             oO = oO + nCB ( 3 )
!           end do !-- iS
!         end do !-- iC
!       end do !-- jC

!       call CO % AllToAll_V ( )

!       oP = 0
!       oI = 0
!       do iG = 1, C % Communicator_3 % Size  /  C % nBricks ( 3 )
!         oC = 0
!         do iB = 1, C % nBricks ( 3 )
!           iR = ( iG - 1 ) * C % nBricks ( 3 )  +  iB  -  1
!           do iPS = 1, C % nSegmentsFrom_3 ( iR )
!             iP = oP + iPS
!             FC % nCoarsen_3 ( iP ) = int ( Incoming ( oI + 1 ) + 0.5_KDR )
!             oI = oI + 1
!             associate ( CP => FC % CoarsenPillar_3 ( iP ) )
!             do iS = 1, CP % nVariables
!               CP % Value ( oC + 1 : oC + nCB ( 3 ), iS )  &
!                 =  Incoming ( oI + 1 : oI + nCB ( 3 ) )
!               oI = oI + nCB ( 3 )
!             end do !-- iS
!             end associate !-- CP
!           end do !-- iPS
!           oC = oC + nCB ( 3 )
!         end do !-- iB
!         oP = oP + C % nSegmentsFrom_3 ( iR )
!       end do !-- iG

!       end associate !-- Outgoing, etc.
!       end associate !-- CO

!     end select !-- iAngular
!     end associate !-- nCB
!     end select !-- C
!     end select !-- PS
!     nullify ( G, Crsn_2, Crsn_3, Vol, SV )

!   end subroutine ComposePillars


!   subroutine CoarsenPillars ( FC, iAngular )

!     class ( FluidCentralTemplate ), intent ( inout ), target :: &
!       FC
!     integer ( KDI ), intent ( in ) :: &
!       iAngular

!     integer ( KDI ) :: &
!       oC, &  !-- oCell
!       iP, &  !-- iPillar
!       iA, &  !-- iAverage
!       iV     !-- iVariable
!     integer ( KDI ), dimension ( : ), pointer :: &
!       nCoarsen
!     real ( KDR ) :: &
!       Volume_A, &
!       Density_A
!     type ( StorageForm ), dimension ( : ), pointer :: &
!       CP

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_C_Template )

!     select case ( iAngular )
!     case ( 2 )
!       if ( .not. C % Communicator_2 % Initialized ) &
!         return
!       nCoarsen => FC % nCoarsen_2
!       CP => FC % CoarsenPillar_2
!     case ( 3 )
!       if ( .not. C % Communicator_3 % Initialized ) &
!         return
!       nCoarsen => FC % nCoarsen_3
!       CP => FC % CoarsenPillar_3
!     end select !-- iAngular

!     do iP = 1, size ( CP )
!       oC = 0
!       associate ( Volume => CP ( iP ) % Value ( :, 1 ) )
!       do iA = 1, CP ( iP ) % nValues / nCoarsen ( iP )
!         Volume_A = sum ( Volume ( oC + 1 : oC + nCoarsen ( iP ) ) )
!         do iV = 2, CP ( iP ) % nVariables
!           associate ( Density => CP ( iP ) % Value ( :, iV ) )
!           Density_A = sum ( Volume ( oC + 1 : oC + nCoarsen ( iP ) ) &
!                             *  Density ( oC + 1 : oC + nCoarsen ( iP ) ) )
!           Density ( oC + 1 : oC + nCoarsen ( iP ) ) &
!             =  Density_A / Volume_A
!           end associate !-- Density
!         end do !-- iV
!         oC = oC + nCoarsen ( iP )
!       end do !-- iA
!       end associate !-- Volume
!     end do !-- iP

!     end select !-- C
!     end select !-- PS
!     nullify ( CP, nCoarsen )

!   end subroutine CoarsenPillars


!   subroutine DecomposePillars ( FC, S, iAngular )

!     class ( FluidCentralTemplate ), intent ( inout ) :: &
!       FC
!     class ( StorageForm ), intent ( inout ) :: &
!       S
!     integer ( KDI ), intent ( in ) :: &
!       iAngular

!     integer ( KDI ) :: &
!       oO, oI, &      !-- oOutgoing, oIncoming
!       oP, &          !-- oPillar
!       oC, &          !-- oCell
!       iC, jC, kC, &  !-- iCell, etc.
!       iS, &          !-- iSelected
!       iG, &          !-- iGroup (of pillars/processes)
!       iB, &          !-- iBrick
!       iR, &          !-- iRank
!       iPS, &         !-- iPillarSegment
!       iP, &          !-- iPillar
!       iMomentum_2
!     real ( KDR ), dimension ( :, :, : ), pointer :: &
!       R, &
!       Crsn_2, Crsn_3, &
!       SV
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_Form )
!     G => PS % Geometry ( )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_C_Template )

!     associate ( nCB => C % nCellsBrick )

!     select case ( iAngular )
!     case ( 2 )

!       if ( .not. C % Communicator_2 % Initialized ) &
!         return

!       associate &
!         ( CO => FC % CO_CoarsenBackward_2 )
!       associate &
!         ( Outgoing => CO % Outgoing % Value, &
!           Incoming => CO % Incoming % Value )

!       oP = 0
!       oO = 0
!       do iG = 1, C % Communicator_2 % Size  /  C % nBricks ( 2 )
!         oC = 0
!         do iB = 1, C % nBricks ( 2 )
!           iR = ( iG - 1 ) * C % nBricks ( 2 )  +  iB  -  1
!           do iPS = 1, C % nSegmentsFrom_2 ( iR )
!             iP = oP + iPS
!             associate ( CP => FC % CoarsenPillar_2 ( iP ) )
!             do iS = 2, CP % nVariables
!               Outgoing ( oO + 1 : oO + nCB ( 2 ) )  &
!                 =  CP % Value ( oC + 1 : oC + nCB ( 2 ), iS )
!               oO = oO + nCB ( 2 )
!             end do !-- iS
!             end associate !-- CP
!           end do !-- iPS
!           oC = oC + nCB ( 2 )
!         end do !-- iB
!         oP = oP + C % nSegmentsFrom_2 ( iR )
!       end do !-- iG

!       call CO % AllToAll_V ( )

!       call C % SetVariablePointer &
!              ( G % Value ( :, G % COARSENING ( 2 ) ), Crsn_2 )
!       call C % SetVariablePointer &
!              ( G % Value ( :, G % CENTER_U ( 1 ) ), R )

!       iMomentum_2 = 3  !-- Fluid

!       oI = 0
!       do kC = 1, nCB ( 3 )
!         do iC = 1, nCB ( 1 )
!           if ( Crsn_2 ( iC, 1, kC ) < 2.0_KDR ) &
!             cycle
!           do iS = 1, S % nVariables
!             call C % SetVariablePointer &
!                    ( S % Value ( :, S % iaSelected ( iS ) ), SV )
!             SV ( iC, 1 : nCB ( 2 ), kC )  &
!               =  Incoming ( oI + 1 : oI + nCB ( 2 ) )
!             oI = oI + nCB ( 2 )
!             if ( iS == iMomentum_2 &
!                  .and. R ( iC, 1, kC )  <  FC % RadiusPolarMomentum ) &
!               SV ( iC, 1 : nCB ( 2 ), kC ) = 0.0_KDR
!           end do !-- iS
!         end do !-- iC
!       end do !-- kC

!       end associate !-- Outgoing, etc.
!       end associate !-- CO

!     case ( 3 )

!       if ( .not. C % Communicator_3 % Initialized ) &
!         return

!       associate &
!         ( CO => FC % CO_CoarsenBackward_3 )
!       associate &
!         ( Outgoing => CO % Outgoing % Value, &
!           Incoming => CO % Incoming % Value )

!       oP = 0
!       oO = 0
!       do iG = 1, C % Communicator_3 % Size  /  C % nBricks ( 3 )
!         oC = 0
!         do iB = 1, C % nBricks ( 3 )
!           iR = ( iG - 1 ) * C % nBricks ( 3 )  +  iB  -  1
!           do iPS = 1, C % nSegmentsFrom_3 ( iR )
!             iP = oP + iPS
!             associate ( CP => FC % CoarsenPillar_3 ( iP ) )
!             do iS = 2, CP % nVariables
!               Outgoing ( oO + 1 : oO + nCB ( 3 ) )  &
!                 =  CP % Value ( oC + 1 : oC + nCB ( 3 ), iS )
!               oO = oO + nCB ( 3 )
!             end do !-- iS
!             end associate !-- CP
!           end do !-- iPS
!           oC = oC + nCB ( 3 )
!         end do !-- iB
!         oP = oP + C % nSegmentsFrom_3 ( iR )
!       end do !-- iG

!       call CO % AllToAll_V ( )

!       call C % SetVariablePointer &
!              ( G % Value ( :, G % COARSENING ( 3 ) ), Crsn_3 )

!       oI = 0
!       do jC = 1, nCB ( 2 )
!         do iC = 1, nCB ( 1 )
!           if ( Crsn_3 ( iC, jC, 1 ) < 2.0_KDR ) &
!             cycle
!           do iS = 1, S % nVariables
!             call C % SetVariablePointer &
!                    ( S % Value ( :, S % iaSelected ( iS ) ), SV )
!             SV ( iC, jC, 1 : nCB ( 3 ) )  &
!               =  Incoming ( oI + 1 : oI + nCB ( 3 ) )
!             oI = oI + nCB ( 3 )
!           end do !-- iS
!         end do !-- iC
!       end do !-- jC

!       end associate !-- Outgoing, etc.
!       end associate !-- CO

!     end select !-- iAngular
!     end associate !-- nCB
!     end select !-- C
!     end select !-- PS
!     nullify ( G, R, Crsn_2, Crsn_3, SV )

!   end subroutine DecomposePillars


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
