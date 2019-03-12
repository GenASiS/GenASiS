module FluidBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ComputeGravity_Command
  use ApplyGravity_F__Command
  use Universe_Template

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: FluidBoxForm
  contains
    procedure, private, pass :: &
      Initialize_FB
    generic, public :: &
      Initialize => Initialize_FB
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_FB
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_FB
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      InitializeStep
  end type FluidBoxForm

contains


  subroutine Initialize_FB &
               ( FB, FluidType, GeometryType, Name, GravitySolverTypeOption, &
                 MinCoordinateOption, MaxCoordinateOption, FinishTimeOption, &
                 CourantFactorOption, UniformAccelerationOption, nCellsOption, &
                 nWriteOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      FluidType, &
      GeometryType, &
      Name
    character ( * ), intent ( in ), optional :: &
      GravitySolverTypeOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      UniformAccelerationOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( FB % Type == '' ) &
      FB % Type = 'a FluidBox'

    call FB % InitializeTemplate ( Name )

    call FB % AllocateIntegrator &
           ( )
    call FB % InitializePositionSpace &
           ( GeometryType, &
             GravitySolverTypeOption = GravitySolverTypeOption, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             UniformAccelerationOption = UniformAccelerationOption, &
             nCellsOption = nCellsOption )
    call FB % InitializeFluid &
          ( FluidType )
    call FB % InitializeStep &
          ( GravitySolverTypeOption )

    select type ( I => FB % Integrator )
    class is ( Integrator_C_PS_Form )
      call I % Initialize &
             ( Name, TimeUnitOption = FB % Units % Time, &
               FinishTimeOption = FinishTimeOption, &
               CourantFactorOption = CourantFactorOption, &
               nWriteOption = nWriteOption )
    end select !-- I

  end subroutine Initialize_FB


  impure elemental subroutine Finalize ( FB )

    type ( FluidBoxForm ), intent ( inout ) :: &
      FB

    call FB % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine AllocateIntegrator_FB ( FB )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB

    allocate ( Integrator_C_PS_Form :: FB % Integrator )

  end subroutine AllocateIntegrator_FB


  subroutine InitializePositionSpace &
               ( FB, GeometryType, GravitySolverTypeOption, &
                 MinCoordinateOption, MaxCoordinateOption, &
                 UniformAccelerationOption, nCellsOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      GeometryType
    character ( * ), intent ( in ), optional :: &
      GravitySolverTypeOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      UniformAccelerationOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption

    integer ( KDI ) :: &
      iD  !-- iDimension

    associate ( I => FB % Integrator )

    allocate ( Atlas_SC_Form :: I % PositionSpace )
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( allocated ( FB % BoundaryConditionsFace ) ) then
      do iD = 1, PS % nDimensions
        call PS % SetBoundaryConditionsFace &
               ( FB % BoundaryConditionsFace ( iD ) % Value, &
                 iDimension = iD )
      end do !-- iD
    end if !-- BoundaryConditions

    call PS % CreateChart &
           ( CoordinateUnitOption = FB % Units % Coordinate_PS, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    call GA % Initialize &
           ( PS, GeometryType, &
             GravitySolverTypeOption = GravitySolverTypeOption, &
             UniformAccelerationOption = UniformAccelerationOption )
    call PS % SetGeometry ( GA )

    end select !-- GA
    end select !-- PS
    end associate !-- I 

  end subroutine InitializePositionSpace


  subroutine InitializeFluid ( FB, FluidType, BaryonMassReferenceOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      FluidType
    real ( KDR ), intent ( in ), optional :: &
      BaryonMassReferenceOption

    select type ( I => FB % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Fluid_ASC_Form :: I % Current_ASC )
    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize &
           ( PS, FluidType, FB % Units, &
             BaryonMassReferenceOption = BaryonMassReferenceOption )
    end select !-- FA

    end select !-- PS
    end select !-- I

  end subroutine InitializeFluid


  subroutine InitializeStep ( FB, Name, GravitySolverTypeOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      Name
    character ( * ), intent ( in ), optional :: &
      GravitySolverTypeOption

    select type ( I => FB % Integrator )
    class is ( Integrator_C_PS_Form )

    allocate ( Step_RK2_C_ASC_Form :: I % Step )
    select type ( S => I % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( I, I % Current_ASC, Name )
    if ( present ( GravitySolverTypeOption ) ) &   
      S % ComputeConstraints % Pointer => ComputeGravity
      S % ApplySources % Pointer => ApplyGravity_F
    end select !-- S

    end select !-- I

  end subroutine InitializeStep


end module FluidBox_Form
