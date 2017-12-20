module FluidCentralCore_Form

  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: FluidCentralCoreForm
    real ( KDR ) :: &
      GravityFactor
    logical ( KDL ) :: &
      Dimensionless = .false.
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetWriteTimeInterval
    procedure, public, pass :: &
      PrepareCycle
    procedure, private, pass :: &
      ComputeTimeStepLocal
    procedure, public, pass :: &
      ComputeTimeStep_G_ASC
  end type FluidCentralCoreForm

    class ( FluidCentralCoreForm ), private, pointer :: &
      FluidCentralCore => null ( )

    private :: &
      ApplySources

      private :: &
        ApplyGravityMomentum, &
        LocalMax, &
        ComputeTimeStep_G_CSL

contains


  subroutine Initialize &
               ( FCC, Name, FluidType, GeometryType, &
                 DimensionlessOption, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, GravityFactorOption, &
                 LimiterParameterOption, nWriteOption )

    class ( FluidCentralCoreForm ), intent ( inout ), target :: &
      FCC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      GravityFactorOption, &
      LimiterParameterOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax, &
      RadialRatio, &
      FinishTime
    type ( MeasuredValueForm ) :: &
      TimeUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit


    if ( FCC % Type == '' ) &
      FCC % Type = 'a FluidCentralCore'

    FluidCentralCore => FCC

    FCC % Dimensionless = .false.
    if ( present ( DimensionlessOption ) ) &
      FCC % Dimensionless = DimensionlessOption


    !-- PositionSpace

    allocate ( Atlas_SC_CC_Form :: FCC % PositionSpace )
    select type ( PS => FCC % PositionSpace )
    class is ( Atlas_SC_CC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( FCC % Dimensionless ) then

      call PS % CreateChart_CC ( )

    else

      RadiusCore   =   40.0_KDR  *  UNIT % KILOMETER
      RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
      RadialRatio  =  6.5_KDR

      CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

      call PS % CreateChart_CC &
             ( CoordinateUnitOption = CoordinateUnit, &
               RadiusCoreOption = RadiusCore, &
               RadiusMaxOption = RadiusMax, &
               RadialRatioOption = RadialRatio )

    end if !-- Dimensionless

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    if ( FCC % Dimensionless ) then
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'MULTIPOLE', &
               GravitationalConstantOption = 1.0_KDR )
    else
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'MULTIPOLE' )
    end if !-- Dimensionless

    call PS % SetGeometry ( )

    FCC % GravityFactor = 0.1_KDR
    if ( present ( GravityFactorOption ) ) &
      FCC % GravityFactor = GravityFactorOption
    call PROGRAM_HEADER % GetParameter &
           ( FCC % GravityFactor, 'GravityFactor' )


    !-- Fluid

    allocate ( Fluid_ASC_Form :: FCC % Current_ASC )
    select type ( FA => FCC % Current_ASC )
    class is ( Fluid_ASC_Form )

    if ( FCC % Dimensionless ) then
      call FA % Initialize &
             ( PS, FluidType, &
               LimiterParameterOption = LimiterParameterOption )
    else
      call FA % Initialize &
             ( PS, FluidType, &
               BaryonMassUnitOption &
                 =  UNIT % ATOMIC_MASS_UNIT, &
               NumberUnitOption &
                 =  UNIT % SOLAR_BARYON_NUMBER, &
               EnergyUnitOption &
                 =  UNIT % SOLAR_MASS  *  UNIT % SPEED_OF_LIGHT **2, &
               MomentumUnitOption &
                 =  UNIT % SOLAR_MASS  *  UNIT % SPEED_OF_LIGHT, &
               AngularMomentumUnitOption &
                 =  UNIT % SOLAR_KERR_PARAMETER, &
               BaryonMassReferenceOption = CONSTANT % ATOMIC_MASS_UNIT, &
               LimiterParameterOption = LimiterParameterOption )
    end if


    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FCC % Step )
    select type ( S => FCC % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources
    end select !-- S


    !-- Template

    if ( trim ( GeometryType ) == 'NEWTONIAN' ) then
      if ( .not. allocated ( FCC % TimeStepLabel ) ) &
        allocate ( FCC % TimeStepLabel ( 2 ) )
      FCC % TimeStepLabel ( 1 )  =  'Fluid'
      FCC % TimeStepLabel ( 2 )  =  'Gravity'
    end if

    if ( .not. FCC % Dimensionless ) &
      TimeUnit = UNIT % SECOND
    if ( present ( TimeUnitOption ) ) &
      TimeUnit = TimeUnitOption

    FinishTime = 1.0_KDR * TimeUnit
    if ( present ( FinishTimeOption ) ) &
      FinishTime = FinishTimeOption

    call FCC % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, &
             FinishTimeOption = FinishTime, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )
    call Show ( FCC % GravityFactor, 'GravityFactor' )
    call Show ( FCC % Dimensionless, 'Dimensionless' )


    !-- Cleanup

    end select !-- FA
    end select !-- GA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCC )

    type ( FluidCentralCoreForm ), intent ( inout ) :: &
      FCC

    call FCC % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetWriteTimeInterval ( I )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
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
    class ( Fluid_D_Form ), pointer :: &
      F

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    G => PS % Geometry ( )

    select type ( Chart => PS % Chart )
    class is ( Chart_SL_Template )

    associate ( C => PS % Communicator ) 

    !-- Find max velocity
    allocate ( CO )
    call CO % Initialize ( C, nOutgoing = [ 1 ], nIncoming = [ C % Size ] )
    CO % Outgoing % Value ( 1 ) &
      =  LocalMax ( Chart % IsProperCell, &
                    abs ( F % Value ( :, F % VELOCITY_U ( 1 ) ) ) ) 
    call CO % Gather ( )
    VelocityMax  =  maxval ( CO % Incoming % Value )
    iProcess     =  maxloc ( CO % Incoming % Value, dim = 1 )  -  1
    deallocate ( CO )

    if ( VelocityMax == 0.0_KDR ) &
      return

    !-- Find radius of max velocity
    allocate ( CO )
    call CO % Initialize &
           ( C, nOutgoing = [ 1 ], nIncoming = [ 1 ], RootOption = iProcess )
    if ( C % Rank == iProcess ) then
      iRadius  =  maxloc ( abs ( F % Value ( :, F % VELOCITY_U ( 1 ) ) ), &
                           dim = 1 )
      CO % Outgoing % Value ( 1 )  =  G % Value ( iRadius, G % CENTER_U ( 1 ) )
    end if
    call CO % Broadcast ( )
    VelocityMaxRadius = CO % Incoming % Value ( 1 )
    deallocate ( CO )

    !-- Compute average density
    select type ( TI => FA % TallyInterior )
    class is ( Tally_F_D_Form )
    DensityAve  =  F % BaryonMassReference &
                   * TI % Value ( TI % BARYON_NUMBER ) &
                   / ( 4.0 / 3.0  *  CONSTANT % PI  *  VelocityMaxRadius ** 3 )

    !-- Time scales
    TimeScaleVelocityMax &
      =  VelocityMaxRadius  /  VelocityMax
    TimeScaleDensityAve &
      =  ( GA % GravitationalConstant  *  DensityAve ) ** ( -0.5_KDR )

    I % WriteTimeInterval  &
      =  min ( TimeScaleDensityAve, TimeScaleVelocityMax )  /  I % nWrite

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( VelocityMax, Chart % CoordinateUnit ( 1 ) / I % TimeUnit, &
                'VelocityMax', I % IGNORABILITY )
    call Show ( VelocityMaxRadius, Chart % CoordinateUnit ( 1 ), &
                'VelocityMaxRadius', I % IGNORABILITY )
    call Show ( DensityAve, UNIT % IDENTITY, 'DensityAve', &
                I % IGNORABILITY )
    call Show ( TimeScaleDensityAve, I % TimeUnit, 'TimeScaleDensityAve', &
                I % IGNORABILITY )
    call Show ( TimeScaleVelocityMax, I % TimeUnit, 'TimeScaleVelocityMax', &
                I % IGNORABILITY )

    !-- Cleanup
    end select !-- TI
    end associate !-- C
    end select !-- Chart
    end select !-- GA
    end select !-- PS
    end select !-- FA
    nullify ( G, F )

  end subroutine SetWriteTimeInterval


  subroutine PrepareCycle ( I )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      I

    ! class ( Fluid_D_Form ), pointer :: &
    !   F

    ! select type ( FA => I % Current_ASC )
    ! class is ( Fluid_ASC_Form )

    ! select type ( PS => I % PositionSpace )
    ! class is ( Atlas_SC_Form )

    ! select type ( GA => PS % Geometry_ASC )
    ! class is ( Geometry_ASC_Form )

    ! F => FA % Fluid_D ( )

    ! call GA % ComputeGravity &
    !       ( I % Current_ASC, &
    !         iBaryonMass = F % BARYON_MASS, &
    !         iBaryonDensity = F % CONSERVED_BARYON_DENSITY )

    ! end select !-- GA
    ! end select !-- PS
    ! end select !-- FA
    ! nullify ( F )

  end subroutine PrepareCycle


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( FluidCentralCoreForm ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    call I % ComputeTimeStep_C_ASC &
           ( TimeStepCandidate ( 1 ), I % Current_ASC )

    call I % ComputeTimeStep_G_ASC &
           ( TimeStepCandidate ( 2 ) )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeTimeStep_G_ASC ( FCC, TimeStepCandidate )

    class ( FluidCentralCoreForm ), intent ( inout ), target :: &
      FCC
    real ( KDR ), intent ( inout ) :: &
      TimeStepCandidate

    class ( Geometry_N_Form ), pointer :: &
      G

    select type ( PS => FCC % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    if ( trim ( GA % GeometryType ) /= 'NEWTONIAN' ) &
      return

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    G  =>  GA % Geometry_N ( )

    call ComputeTimeStep_G_CSL &
           ( CSL % IsProperCell, &
             G % Value ( :, G % GRAVITATIONAL_ACCELERATION_D ( 1 ) ), &
             G % Value ( :, G % GRAVITATIONAL_ACCELERATION_D ( 2 ) ), &
             G % Value ( :, G % GRAVITATIONAL_ACCELERATION_D ( 3 ) ), &
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
    nullify ( G )

  end subroutine ComputeTimeStep_G_ASC


  subroutine ApplySources &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iMomentum
    type ( TimerForm ), pointer :: &
      Timer
    class ( Geometry_N_Form ), pointer :: &
      G

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    call ApplyCurvilinear_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    select type ( PS => FluidCentralCore % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    G  =>  GA % Geometry_N ( )

    select type ( F => Fluid )
    class is ( Fluid_D_Form )

    call GA % ComputeGravity &
           ( FluidCentralCore % Current_ASC, &
             iBaryonMass = F % BARYON_MASS, &
             iBaryonDensity = F % CONSERVED_BARYON_DENSITY )

    select type ( FS => Sources_F )
    class is ( Sources_F_Form )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )
      do iD = 1, C % nDimensions
        call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( iD ), &
                      iMomentum )
        call ApplyGravityMomentum &
               ( F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
                 G % Value ( :, G % GRAVITATIONAL_ACCELERATION_D ( iD ) ), &
                 TimeStep, S % B ( iStage ), &
                 Increment % Value ( :, iMomentum ), & 
                 FS % Value ( :, FS % GRAVITATIONAL_S_D ( iD ) ) )
      end do !-- iD
    end select !-- C

    end select !-- FS
    end select !-- F
    end select !-- GA
    end select !-- PS
    nullify ( G )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplySources


  subroutine ApplyGravityMomentum ( M, N, GradPhi, dt, Weight_RK, K, S )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      GradPhi
    real ( KDR ) :: &
      dt, &
      Weight_RK
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      K, &
      S

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues

    nValues = size ( K )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      K ( iV )  =  K ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                *  dt
      S ( iV )  =  S ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                *  Weight_RK
    end do
    !$OMP end parallel do

  end subroutine ApplyGravityMomentum


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
