module FluidCentralCore_Form

  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: FluidCentralCoreForm
    type ( Storage_ASC_Form ), allocatable :: &
      PoissonStorage
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, public, nopass :: &
      ComputeGravity
  end type FluidCentralCoreForm

contains


  subroutine Initialize &
               ( FCC, Name, FluidType, GeometryType, &
                 GravitySolverTypeOption, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
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
    class ( Geometry_N_Form ), pointer :: &
      G_N

    if ( FCC % Type == '' ) &
      FCC % Type = 'a FluidCentralCore'

    !-- PositionSpace

    allocate ( Atlas_SC_CC_Form :: FCC % PositionSpace )
    select type ( PS => FCC % PositionSpace )
    class is ( Atlas_SC_CC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )
    call PS % CreateChart_CC ( )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    call GA % Initialize ( PS, GeometryType )

    call PS % SetGeometry ( )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: FCC % Current_ASC )
    select type ( FA => FCC % Current_ASC )  !-- FluidAtlas
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, FluidType )

    !-- Gravity

    if ( trim ( GeometryType ) == 'NEWTONIAN' ) then
      if ( present ( GravitySolverTypeOption ) ) then

        MaxDegree = 10
        call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

        G_N => GA % Geometry_N ( )
        allocate ( FCC % PoissonStorage )
        associate ( PSA => FCC % PoissonStorage )
        call PSA % Initialize &
               ( GA, NameShort = 'PoissonStorage', &
                 iaSelectedOption = [ G_N % GRAVITATIONAL_POTENTIAL ] )
        end associate !-- PSA

        allocate ( FCC % Poisson )
        associate ( P => FCC % Poisson )
        call P % Initialize &
               ( PS, SolverType = GravitySolverTypeOption, &
                 MaxDegreeOption = MaxDegree, &
                 nEquationsOption = 1 )
        end associate !-- P

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
           ( Name, FinishTimeOption = FinishTimeOption )

    !-- Cleanup

    end select !-- FA
    end select !-- GA
    end select !-- PS
    nullify ( G_N )

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCC )

    type ( FluidCentralCoreForm ), intent ( inout ) :: &
      FCC

    if ( allocated ( FCC % Poisson ) ) &
      deallocate ( FCC % Poisson )
    if ( allocated ( FCC % PoissonStorage ) ) &
      deallocate ( FCC % PoissonStorage )

    call FCC % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine ComputeGravity ( )

!    call P % Solve ( PFT % Solution, PFT % Source )

  end subroutine ComputeGravity


end module FluidCentralCore_Form
