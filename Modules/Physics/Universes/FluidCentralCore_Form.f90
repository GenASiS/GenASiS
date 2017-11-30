module FluidCentralCore_Form

  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: FluidCentralCoreForm
    type ( Storage_ASC_Form ), allocatable :: &
      Storage_ASC
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, public, nopass :: &
      ComputeGravity
  end type FluidCentralCoreForm

    class ( FluidCentralCoreForm ), private, pointer :: &
      FluidCentralCore => null ( )

    private :: &
      ComputeGravitySource

contains


  subroutine Initialize &
               ( FCC, Name, FluidType, GeometryType, &
                 GravitySolverTypeOption, TimeUnitOption, &
                 FinishTimeOption, CourantFactorOption, nWriteOption )

    class ( FluidCentralCoreForm ), intent ( inout ), target :: &
      FCC
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    character ( * ), intent ( in ), optional :: &
      GravitySolverTypeOption
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
      RadialRatio
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    class ( Geometry_N_Form ), pointer :: &
      G_N

    if ( FCC % Type == '' ) &
      FCC % Type = 'a FluidCentralCore'

    FluidCentralCore => FCC


    !-- PositionSpace

    allocate ( Atlas_SC_CC_Form :: FCC % PositionSpace )
    select type ( PS => FCC % PositionSpace )
    class is ( Atlas_SC_CC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    RadiusCore   =   40.0_KDR  *  UNIT % KILOMETER
    RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
    RadialRatio  =  6.5_KDR

    CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

    call PS % CreateChart_CC &
           ( CoordinateUnitOption = CoordinateUnit, &
             RadiusCoreOption = RadiusCore, &
             RadiusMaxOption = RadiusMax, &
             RadialRatioOption = RadialRatio )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    call GA % Initialize ( PS, GeometryType )

    call PS % SetGeometry ( )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: FCC % Current_ASC )
    select type ( FA => FCC % Current_ASC )  !-- FluidAtlas
    class is ( Fluid_ASC_Form )

    select case ( trim ( GeometryType ) )
    case ( 'GALILEAN' )
      call FA % Initialize ( PS, FluidType )
    case default
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
    end select

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
        call Show ( 'NEWTONIAN geometry required GravitySolverType', &
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
    end select !-- S

    !-- Template

    call FCC % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption = FinishTimeOption, &
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


  subroutine ComputeGravity ( )

    type ( VariableGroupForm ), pointer :: &
      S
    class ( Fluid_D_Form ), pointer :: &
      F

    associate &
      ( PA  => FluidCentralCore % Poisson_ASC, &
        SA => FluidCentralCore % Storage_ASC )

    select type ( FA => FluidCentralCore % Current_ASC )  !-- FluidAtlas
    class is ( Fluid_ASC_Form )

    S  =>  SA % Storage ( )
    F  =>  FA % Fluid_D ( )

    call ComputeGravitySource &
           ( F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
             S % Value ( :, S % iaSelected ( 1 ) ) )

    call PA % Solve ( SA, SA )

    end select !-- FA
    end associate !-- PA, etc.
    nullify ( S, F )

  end subroutine ComputeGravity


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

    FourPi_G  =  4.0_KDR  *  CONSTANT % PI  *  CONSTANT % GRAVITATIONAL

    nValues = size ( S )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      S ( iV )  =  FourPi_G  *  M ( iV )  *  N ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeGravitySource


end module FluidCentralCore_Form
