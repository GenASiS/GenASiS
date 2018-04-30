module FluidBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ApplyGravity_F__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: FluidBoxForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FluidBoxForm

contains


  subroutine Initialize &
               ( FB, Name, FluidType, GeometryType, &
                 BoundaryConditionsFaceOption, GravitySolverTypeOption, &
                 MinCoordinateOption, MaxCoordinateOption, TimeUnitOption, &
                 FinishTimeOption, CourantFactorOption, &
                 UniformAccelerationOption, nCellsOption, nWriteOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    type ( Character_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      BoundaryConditionsFaceOption
    character ( * ), intent ( in ), optional :: &
      GravitySolverTypeOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      UniformAccelerationOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iD  !-- iDimension

    if ( FB % Type == '' ) &
      FB % Type = 'a FluidBox'

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: FB % PositionSpace )
    select type ( PS => FB % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( present ( BoundaryConditionsFaceOption ) ) then
      do iD = 1, PS % nDimensions
        call PS % SetBoundaryConditionsFace &
               ( BoundaryConditionsFaceOption ( iD ) % Value, &
                 iDimension = iD )
      end do !-- iD
    end if

    call PS % CreateChart &
           ( MinCoordinateOption = MinCoordinateOption, &
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

    !-- Fluid

    allocate ( Fluid_ASC_Form :: FB % Current_ASC )
    select type ( FA => FB % Current_ASC )  !-- FluidAtlas
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, FluidType )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: FB % Step )
    select type ( S => FB % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( FA, Name )
    if ( present ( GravitySolverTypeOption ) ) &   
      S % ApplySources % Pointer => ApplyGravity_F
    end select !-- S

    !-- Template

    call FB % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption = FinishTimeOption )

    !-- Cleanup

    end select !-- FA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( FB )

    type ( FluidBoxForm ), intent ( inout ) :: &
      FB

    call FB % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


end module FluidBox_Form
