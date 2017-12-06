module FluidCentralCore_Form

  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: FluidCentralCoreForm
    real ( KDR ) :: &
      GravitationalConstant
    logical ( KDL ) :: &
      Dimensionless = .false.
    type ( Storage_ASC_Form ), allocatable :: &
      Storage_ASC
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &  !-- 3
      SetWriteTimeInterval
    procedure, public, pass :: &
      PrepareCycle
  end type FluidCentralCoreForm

    class ( FluidCentralCoreForm ), private, pointer :: &
      FluidCentralCore => null ( )

    private :: &
      ApplySources, &
      ComputeGravity

      private :: &
        ApplyGravityMomentum, &
        ComputeGravitySource, &
        ComputeGravitationalForce, &
        LocalMax

contains


  subroutine Initialize &
               ( FCC, Name, FluidType, GeometryType, &
                 GravitySolverTypeOption, DimensionlessOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 nWriteOption )

    class ( FluidCentralCoreForm ), intent ( inout ), target :: &
      FCC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    character ( * ), intent ( in ), optional :: &
      GravitySolverTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      MaxDegree
    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax, &
      RadialRatio, &
      FinishTime
    type ( MeasuredValueForm ) :: &
      TimeUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    class ( Geometry_N_Form ), pointer :: &
      G_N

    if ( FCC % Type == '' ) &
      FCC % Type = 'a FluidCentralCore'

    FluidCentralCore => FCC

    if ( present ( DimensionlessOption ) ) &
      FCC % Dimensionless = DimensionlessOption

    if ( FCC % Dimensionless ) then
      FCC % GravitationalConstant  =  1.0_KDR
    else
      FCC % GravitationalConstant  =  CONSTANT % GRAVITATIONAL
    end if


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
    call GA % Initialize ( PS, GeometryType )

    call PS % SetGeometry ( )


    !-- Fluid

    allocate ( Fluid_ASC_Form :: FCC % Current_ASC )
    select type ( FA => FCC % Current_ASC )
    class is ( Fluid_ASC_Form )

    if ( FCC % Dimensionless ) then
      call FA % Initialize ( PS, FluidType )
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
               BaryonMassReferenceOption = CONSTANT % ATOMIC_MASS_UNIT )
    end if


    !-- Gravity

    if ( trim ( GeometryType ) == 'NEWTONIAN' ) then
      if ( present ( GravitySolverTypeOption ) ) then

        MaxDegree = 10
        call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

        G_N => GA % Geometry_N ( )
        allocate ( FCC % Storage_ASC )
        associate ( SA => FCC % Storage_ASC )
        call SA % Initialize &
               ( GA, NameShort = 'PoissonStorage', &
                 iaSelectedOption = [ G_N % GRAVITATIONAL_POTENTIAL ] )
        end associate !-- SA

        allocate ( FCC % Poisson_ASC )
        associate ( PA => FCC % Poisson_ASC )
        call PA % Initialize &
               ( PS, SolverType = GravitySolverTypeOption, &
                 MaxDegreeOption = MaxDegree, &
                 nEquationsOption = 1 )
        end associate !-- PA

        select type ( TI => FA % TallyInterior )
        class is ( Tally_F_D_Form )
          TI % ComputeGravitationalPotential => ComputeGravity
        end select !-- TI

      else
        call Show ( 'NEWTONIAN geometry requires GravitySolverType', &
                    CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call Show ( 'FluidCentralCore_Form', 'module', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
    end if


    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FCC % Step )
    select type ( S => FCC % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    S % ApplySources % Pointer => ApplySources
    end select !-- S


    !-- Template

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


    !-- Cleanup

    end select !-- FA
    end select !-- GA
    end select !-- PS
    nullify ( G_N )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCC )

    type ( FluidCentralCoreForm ), intent ( inout ) :: &
      FCC

    if ( allocated ( FCC % Poisson_ASC ) ) &
      deallocate ( FCC % Poisson_ASC )
    if ( allocated ( FCC % Storage_ASC ) ) &
      deallocate ( FCC % Storage_ASC )

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
      CO % Outgoing % Value ( 1 )  =  G % Value ( iRadius, G % CENTER ( 1 ) )
    end if
    call CO % Broadcast ( )
    VelocityMaxRadius = CO % Incoming % Value ( 1 )
    deallocate ( CO )

    !-- Compute average density
    select type ( TI => FA % TallyInterior )
    class is ( Tally_F_D_Form )
    DensityAve  =  F % BaryonMassReference &
                   * TI % Value ( TI % BARYON_NUMBER ) &
                   / VelocityMaxRadius ** 3

    !-- Time scales
    TimeScaleVelocityMax &
      =  VelocityMaxRadius  /  VelocityMax
    TimeScaleDensityAve &
      =  ( I % GravitationalConstant  *  DensityAve ) ** ( -0.5_KDR )

    I % WriteTimeInterval  &
      =  min ( TimeScaleDensityAve, TimeScaleVelocityMax )  /  I % nWrite

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( VelocityMax, Chart % CoordinateUnit ( 1 ) / I % TimeUnit, &
                'VelocityMax', I % IGNORABILITY )
    call Show ( VelocityMaxRadius, Chart % CoordinateUnit ( 1 ), &
                'VelocityMaxRadius', I % IGNORABILITY )
    call Show ( DensityAve, UNIT % MASS_DENSITY_CGS, 'DensityAve', &
                I % IGNORABILITY )
    call Show ( TimeScaleDensityAve, I % TimeUnit, 'TimeScaleDensityAve', &
                I % IGNORABILITY )
    call Show ( TimeScaleVelocityMax, I % TimeUnit, 'TimeScaleVelocityMax', &
                I % IGNORABILITY )

    !-- Cleanup
    end select !-- TI
    end associate !-- C
    end select !-- Chart
    end select !-- PS
    end select !-- FA
    nullify ( G, F )

  end subroutine SetWriteTimeInterval


  subroutine PrepareCycle ( I )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      I

    call ComputeGravity ( )

  end subroutine PrepareCycle


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

    select type ( FS => Sources_F )
    class is ( Sources_F_Form )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )
      do iD = 1, C % nDimensions
        call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( iD ), &
                      iMomentum )
        call ApplyGravityMomentum &
               ( G % Value ( :, G % GRAVITATIONAL_FORCE_D ( iD ) ), &
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


  subroutine ComputeGravity ( )

    type ( VariableGroupForm ), pointer :: &
      S
    class ( Geometry_N_Form ), pointer :: &
      G
    class ( Fluid_D_Form ), pointer :: &
      F

    integer ( KDI ) :: &
      iD  !-- iDimension

    associate &
      ( PA  => FluidCentralCore % Poisson_ASC, &
        SA => FluidCentralCore % Storage_ASC )

    select type ( FA => FluidCentralCore % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => FluidCentralCore % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    S  =>  SA % Storage ( )
    G  =>  GA % Geometry_N ( )
    F  =>  FA % Fluid_D ( )

    call ComputeGravitySource &
           ( F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
             S % Value ( :, S % iaSelected ( 1 ) ) )

    call PA % Solve ( SA, SA )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )
      do iD = 1, C % nDimensions
        call GA % Gradient % Compute ( C, S, iDimension = iD )
        call ComputeGravitationalForce & 
               ( F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
                 GA % Gradient % Output % Value ( :, 1 ), &
                 G % Value ( :, G % GRAVITATIONAL_FORCE_D ( iD ) ) ) 
      end do !-- iD
    end select !-- C

    end select !-- GA
    end select !-- PS
    end select !-- FA
    end associate !-- PA, etc.
    nullify ( S, G, F )

  end subroutine ComputeGravity


  subroutine ApplyGravityMomentum ( F, dt, Weight_RK, K, S )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      F
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
      K ( iV )  =  K ( iV )  +  F ( iV )  *  dt
      S ( iV )  =  S ( iV )  +  F ( iV )  *  Weight_RK
    end do
    !$OMP end parallel do

  end subroutine ApplyGravityMomentum


  subroutine ComputeGravitySource ( M, N, S )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      S

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      FourPi_G

    FourPi_G  =  4.0_KDR  *  CONSTANT % PI  &
                 *  FluidCentralCore % GravitationalConstant

    nValues = size ( S )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      S ( iV )  =  FourPi_G  *  M ( iV )  *  N ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeGravitySource


  subroutine ComputeGravitationalForce ( M, N, GradPhi, F )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      GradPhi
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      F

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues

    nValues = size ( F )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      F ( iV )  =  - M ( iV )  *  N ( iV )  *  GradPhi ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeGravitationalForce


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


end module FluidCentralCore_Form
