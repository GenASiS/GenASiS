module FluidCentralCore_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use FluidCentral_Template

  implicit none
  private

  type, public, extends ( FluidCentralTemplate ) :: FluidCentralCoreForm
    real ( KDR ) :: &
      GravityFactor
  contains
    procedure, private, pass :: &
      Initialize_FCC
    generic, public :: &
      Initialize => Initialize_FCC
    final :: &
      Finalize
    procedure, private, pass :: &
      InitializeAtlas
    procedure, private, pass :: &
      InitializeGeometry
    procedure, public, pass :: &
      SetCentralCore
    procedure, public, pass :: &
      SetCoarsening
    procedure, public, nopass :: &
      CoarsenSingularities
    procedure, public, pass :: &
      ComputeTimeStep_G_ASC
  end type FluidCentralCoreForm

      private :: &
        SetCheckpointTimeInterval, &
        ComputeTimeStepLocal, &
        LocalMax, &
        ComputeTimeStep_G_CSL

contains


  subroutine Initialize_FCC &
               ( FCC, FluidType, GeometryType, Name, DimensionlessOption, &
                 FluidUseDeviceOption, GeometryUseDeviceOption, &
                 FinishTimeOption, CourantFactorOption, GravityFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, &
                 RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
                 nCellsPolarOption, nWriteOption )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      FCC
    character ( * ), intent ( in )  :: &
      FluidType, &
      GeometryType, &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption, &
      FluidUseDeviceOption, &
      GeometryUseDeviceOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      GravityFactorOption, &
      LimiterParameterOption, &
      ShockThresholdOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( FCC % Type == '' ) &
      FCC % Type = 'a FluidCentralCore'

    if ( any ( trim ( GeometryType ) &
                 == [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
    then
      if ( .not. allocated ( FCC % TimeStepLabel ) ) &
        allocate ( FCC % TimeStepLabel ( 2 ) )
      FCC % TimeStepLabel ( 1 )  =  'Fluid_Advection'
      FCC % TimeStepLabel ( 2 )  =  'Gravity'
    end if

    call FCC % InitializeTemplate_FC &
           ( FluidType, GeometryType, Name, &
             DimensionlessOption = DimensionlessOption, &
             FluidUseDeviceOption = FluidUseDeviceOption, &
             GeometryUseDeviceOption = GeometryUseDeviceOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             LimiterParameterOption = LimiterParameterOption, &
             ShockThresholdOption = ShockThresholdOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nWriteOption = nWriteOption, &
             nCellsPolarOption = nCellsPolarOption )

    call FCC % SetCentralCore ( GeometryType, GravityFactorOption )

  end subroutine Initialize_FCC


  subroutine InitializeAtlas &
               ( FC, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      FC
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusCore, &
      RadialRatio

    associate ( I => FC % Integrator )

    allocate ( Atlas_SC_CC_Form :: I % PositionSpace )
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_CC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( FC % Dimensionless ) then

      call PS % CreateChart_CC ( )

    else

      RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
      RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
      RadialRatio  =  2.4_KDR

      call PS % CreateChart_CC &
             ( CoordinateUnitOption = FC % Units % Coordinate_PS, &
               RadiusCoreOption = RadiusCore, &
               RadiusMaxOption = RadiusMax, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

      FC % RadiusPolarMomentum  =  8.0_KDR  *  UNIT % KILOMETER
      call PROGRAM_HEADER % GetParameter &
             ( FC % RadiusPolarMomentum, 'RadiusPolarMomentum' )

    end if !-- Dimensionless

    end select !-- PS
    end associate !-- I

  end subroutine InitializeAtlas


  subroutine InitializeGeometry &
               ( FC, GA, PS, GeometryType, UsePinnedMemoryOption, &
                 CentralMassOption )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      FC
    type ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      PS
    character ( * ), intent ( in )  :: &
      GeometryType
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    real ( KDR ), intent ( in ), optional :: &
      CentralMassOption

    character ( LDL ) :: &
      GravitySolverType

    GravitySolverType  =  'MULTIPOLE'
    !-- FIXME: XL 16.1.1-5 does not work without association.
    associate ( PH => PROGRAM_HEADER )
!    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )
    call PH % GetParameter ( GravitySolverType, 'GravitySolverType' )
    end associate !-- PH    

    if ( FC % Dimensionless ) then
      call GA % Initialize &
             ( PS, GeometryType, & 
               UsePinnedMemoryOption = UsePinnedMemoryOption, &
               GravitySolverTypeOption = GravitySolverType, &
               GravitationalConstantOption = 1.0_KDR )
    else
      call GA % Initialize &
             ( PS, GeometryType, &
               UsePinnedMemoryOption = UsePinnedMemoryOption, &
               GravitySolverTypeOption = GravitySolverType )
    end if !-- Dimensionless

  end subroutine InitializeGeometry


  subroutine SetCentralCore ( FCC, GeometryType, GravityFactorOption )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      FCC
    character ( * ), intent ( in ) :: &
      GeometryType
    real ( KDR ), intent ( in ), optional :: &
      GravityFactorOption

    select type ( I => FCC % Integrator )
    class is ( Integrator_C_PS_Form )

    I % SetCheckpointTimeInterval => SetCheckpointTimeInterval
    I % ComputeTimeStepLocal => ComputeTimeStepLocal

    select type ( S => I % Step )
    class is ( Step_RK_C_ASC_Template )
      S % CoarsenSingularities => CoarsenSingularities
    end select !-- S

    end select !-- I

    call Show ( 'FluidCentralCore parameter' )
    if ( any ( trim ( GeometryType ) &
                 == [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
    then
      FCC % GravityFactor = 0.7_KDR
      if ( present ( GravityFactorOption ) ) &
        FCC % GravityFactor = GravityFactorOption
      call PROGRAM_HEADER % GetParameter &
             ( FCC % GravityFactor, 'GravityFactor' )
      call Show ( FCC % GravityFactor, 'GravityFactor', FCC % IGNORABILITY )
    end if

  end subroutine SetCentralCore


  impure elemental subroutine Finalize ( FCC )

    type ( FluidCentralCoreForm ), intent ( inout ) :: &
      FCC

    call FCC % FinalizeTemplate_FC ( )

  end subroutine Finalize


  subroutine SetCoarsening ( FC )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      FC

    select type ( I => FC % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_CC_Form )

    !-- Polar coarsening

    if ( PS % nDimensions < 2 ) &
      return
    call Show ( 'SetCoarsening_2', FC % IGNORABILITY )
    call Show ( FC % Name, 'Universe', FC % IGNORABILITY )
    call FC % SetCoarseningTemplate ( iAngular = 2 )

    !-- Azimuthal coarsening

    if ( PS % nDimensions < 3 ) &
      return
    call Show ( 'SetCoarsening_3', FC % IGNORABILITY )
    call Show ( FC % Name, 'Universe', FC % IGNORABILITY )
    call FC % SetCoarseningTemplate ( iAngular = 3 )

    end select !-- PS
    end select !-- I

  end subroutine SetCoarsening


  subroutine CoarsenSingularities ( S, Increment )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( StorageForm ), intent ( inout ) :: &
      Increment

    select type ( FCC => S % Integrator % Universe )
    class is ( FluidCentralCoreForm )

    select type ( I => S % Integrator )
    class is ( IntegratorTemplate )

    select type ( PS => FCC % Integrator % PositionSpace )
    class is ( Atlas_SC_CC_Form )

    select case ( mod ( FCC % Integrator % iCycle, 2 ) )
    case ( 0 )
      if ( PS % nDimensions > 2 ) &
        call FCC % CoarsenSingularityTemplate ( Increment, iAngular = 3 )
      if ( PS % nDimensions > 1 ) &
        call FCC % CoarsenSingularityTemplate ( Increment, iAngular = 2 )
    case ( 1 )
      if ( PS % nDimensions > 1 ) &
        call FCC % CoarsenSingularityTemplate ( Increment, iAngular = 2 )
      if ( PS % nDimensions > 2 ) &
        call FCC % CoarsenSingularityTemplate ( Increment, iAngular = 3 )
    end select
  
    end select !-- PS
    end select !-- I
    end select !-- FCC

  end subroutine CoarsenSingularities


  subroutine SetCheckpointTimeInterval ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iRadius
    real ( KDR ) :: &
      GravitationalConstant, &
      BaryonMass, &
      VelocityMax, &
      VelocityMaxRadius, &
      NumberDensityMax, &
      TimeScaleVelocityMax, &
      TimeScaleDensityMax
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G_SA
    class ( Fluid_D_Form ), pointer :: &
      F_SA

    select type ( FC => I % Universe )
    class is ( FluidCentralCoreForm )

    GravitationalConstant  =  CONSTANT % GRAVITATIONAL
               BaryonMass  =  CONSTANT % ATOMIC_MASS_UNIT 
    if ( FC % Dimensionless ) then
      GravitationalConstant  =  1.0_KDR
                 BaryonMass  =  1.0_KDR
    end if

    G_SA  =>  FC % PositionSpace_SA % Geometry ( )
    F_SA  =>  FC % Fluid_ASC_SA % Fluid_D ( )

    iRadius &
      =  maxloc ( abs ( F_SA % Value ( :,  F_SA % VELOCITY_U ( 1 ) ) ), &
                  dim = 1 )
    VelocityMaxRadius &
      =  G_SA % Value ( iRadius, G_SA % CENTER_U ( 1 ) )
    VelocityMax  &
      =  max ( maxval ( abs ( F_SA % Value ( :, F_SA % VELOCITY_U ( 1 ) ) ) ), &
               sqrt ( tiny ( 0.0_KDR ) ) )  
    NumberDensityMax  &
      =  max ( maxval ( F_SA % Value ( :, F_SA % COMOVING_BARYON_DENSITY ) ), &
               sqrt ( tiny ( 0.0_KDR ) ) )

    TimeScaleVelocityMax &
      =  VelocityMaxRadius  /  VelocityMax
    TimeScaleDensityMax &
      =  ( GravitationalConstant * BaryonMass * NumberDensityMax ) &
         ** ( -0.5_KDR )

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( VelocityMaxRadius, FC % Units % Coordinate_PS ( 1 ), &
                'VelocityMaxRadius', I % IGNORABILITY )
    call Show ( VelocityMax, FC % Units % Velocity_U ( 1 ), &
                'VelocityMax', I % IGNORABILITY )
    call Show ( NumberDensityMax, FC % Units % NumberDensity, &
                'NumberDensityMax', I % IGNORABILITY )
    call Show ( BaryonMass * NumberDensityMax, FC % Units % MassDensity, &
                'MassDensityMax', I % IGNORABILITY )
    call Show ( TimeScaleVelocityMax, I % TimeUnit, &
                'TimeScaleVelocityMax', I % IGNORABILITY )
    call Show ( TimeScaleDensityMax, I % TimeUnit, &
                'TimeScaleDensityMax', I % IGNORABILITY )

    allocate ( CO )
    associate ( C => I % PositionSpace % Communicator ) 
    call CO % Initialize &
           ( C, nOutgoing = [ 1 ], nIncoming = [ 1 ], &
             RootOption = CONSOLE % DisplayRank )
    end associate !-- C

    CO % Outgoing % Value ( 1 )  &
      =  min ( TimeScaleVelocityMax, TimeScaleDensityMax )  /  I % nWrite

    call CO % Broadcast ( )

    I % CheckpointTimeInterval  =  CO % Incoming % Value ( 1 )

    end select !-- FC

    nullify ( G_SA, F_SA )

  end subroutine SetCheckpointTimeInterval


  subroutine ComputeTimeStep_G_ASC ( FCC, TimeStepCandidate )

    class ( FluidCentralCoreForm ), intent ( inout ), target :: &
      FCC
    real ( KDR ), intent ( inout ) :: &
      TimeStepCandidate

    class ( Geometry_N_Form ), pointer :: &
      G

    associate ( I => FCC % Integrator )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    if ( all ( trim ( GA % GeometryType ) &
                 /= [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
      return

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    G  =>  GA % Geometry_N ( )

    call ComputeTimeStep_G_CSL &
           ( CSL % IsProperCell, &
             G % Value ( :, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
             G % Value ( :, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
             G % Value ( :, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
             G % Value ( :, G % METRIC_UU_22 ), &
             G % Value ( :, G % METRIC_UU_33 ), &
             G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
             G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), & 
             G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
             G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), & 
             G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
             CSL % nDimensions, TimeStepCandidate )

    TimeStepCandidate  =  FCC % GravityFactor * TimeStepCandidate

    end select !-- CSL
    end select !-- GA
    end select !-- PS
    end associate !-- I
    nullify ( G )

  end subroutine ComputeTimeStep_G_ASC


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( IntegratorTemplate ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    select type ( I )
    class is ( Integrator_C_PS_Form )
      call I % ComputeTimeStep_C_ASC &
             ( TimeStepCandidate ( 1 ), I % Current_ASC )
    end select !-- I

    select type ( FCC => I % Universe )
    class is ( FluidCentralCoreForm )
      call FCC % ComputeTimeStep_G_ASC &
             ( TimeStepCandidate ( 2 ) )
    end select !-- I

  end subroutine ComputeTimeStepLocal


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


  subroutine ComputeTimeStep_G_CSL &
               ( IsProperCell, GradPhi_1, GradPhi_2, GradPhi_3, &
                 M_UU_22, M_UU_33, dXL_1, dXL_2, dXL_3, dXR_1, dXR_2, dXR_3, &
                 nDimensions, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      GradPhi_1, GradPhi_2, GradPhi_3, &
      M_UU_22, M_UU_33, &
      dXL_1, dXL_2, dXL_3, &
      dXR_1, dXR_2, dXR_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      TimeStepInverse

    nV = size ( dXL_1 )

    select case ( nDimensions )
    case ( 1 )
      TimeStepInverse &
        = maxval ( sqrt ( abs ( GradPhi_1 ) / ( dXL_1 + dXR_1 ) ), &
                   mask = IsProperCell )
    case ( 2 )
      TimeStepInverse &
        = maxval ( sqrt (    abs ( GradPhi_1 ) &
                             / ( dXL_1 + dXR_1 ) &
                          +  abs ( M_UU_22 * GradPhi_2 ) &
                             / ( dXL_2 + dXR_2 ) ), &
                   mask = IsProperCell )
    case ( 3 )
      ! TimeStepInverse &
      !   = maxval ( sqrt (   abs ( GradPhi_1 ) / dX_1 &
      !                     + abs ( M_UU_22 * GradPhi_2 ) / dX_2 &
      !                     + abs ( M_UU_33 * GradPhi_3 ) / dX_3 ) ), &
      !              mask = IsProperCell )
      TimeStepInverse = - huge ( 0.0_KDR )
      !$OMP parallel do private ( iV ) &
      !$OMP reduction ( max : TimeStepInverse )
      do iV = 1, nV
        if ( IsProperCell ( iV ) ) &
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt (   abs ( GradPhi_1 ( iV ) ) &
                                   / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                           + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                                   / ( dXL_2 ( iV ) + dXR_2 ( iV ) ) &
                           + abs ( M_UU_33 ( iV ) * GradPhi_3 ( iV ) ) &
                                   / ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
      end do
      !$OMP end parallel do
    end select !-- nDimensions

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_G_CSL


end module FluidCentralCore_Form
