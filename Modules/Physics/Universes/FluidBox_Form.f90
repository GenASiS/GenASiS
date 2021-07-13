module FluidBox_Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: FluidBoxForm
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
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
      InitializeGravitation
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      InitializeStep
   end type FluidBoxForm


contains


  subroutine Initialize_FB &
               ( FB, FluidType, GravitationType, NameOption, &
                 MinCoordinateOption, MaxCoordinateOption, FinishTimeOption, &
                 nCellsOption, nWriteOption )
!                 CourantFactorOption, UniformAccelerationOption, &

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption!, &
!       CourantFactorOption, &
!       UniformAccelerationOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    character ( LDL ) :: &
      Name

    if ( FB % Type  ==  '' ) &
      FB % Type  =  'a FluidBox'

    Name  =  'FluidBox'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call FB % Universe_H_Form % Initialize ( NameOption = Name )

    call FB % AllocateIntegrator &
           ( )
    call FB % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )
    call FB % InitializeGravitation &
           ( GravitationType )
!              UniformAccelerationOption = UniformAccelerationOption, &
    call FB % InitializeFluid &
           ( FluidType )
    call FB % InitializeStep &
           ( )

    select type ( I  =>  FB % Integrator )
      class is ( Integrator_CS_Form )
    call I % Initialize &
           ( Unit_T_Option = FB % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
!             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )
    end select !-- I

  end subroutine Initialize_FB


  impure elemental subroutine Finalize ( FB )

    type ( FluidBoxForm ), intent ( inout ) :: &
      FB

    if ( allocated ( FB % Units_F ) ) &
      deallocate ( FB % Units_F )

  end subroutine Finalize


  subroutine AllocateIntegrator_FB ( FB )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB

    allocate ( Integrator_CS_Form :: FB % Integrator )

  end subroutine AllocateIntegrator_FB


  subroutine InitializePositionSpace &
               ( FB, MinCoordinateOption, MaxCoordinateOption, nCellsOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption

!     integer ( KDI ) :: &
!       iD  !-- iDimension
    
    associate ( I  =>  FB % Integrator )

    allocate ( Atlas_SCG_Form  ::  I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_Form )

    allocate ( FB % Units_F ( 1 ) )

    call PS % Initialize &
           ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
             NameOption = 'PositionSpace', &
             CoordinateUnitOption = FB % Units_F ( 1 ) % Coordinate_PS, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )

!     if ( allocated ( FB % BoundaryConditionsFace ) ) then
!       do iD = 1, PS % nDimensions
!         call PS % SetBoundaryConditionsFace &
!                ( FB % BoundaryConditionsFace ( iD ) % Value, &
!                  iDimension = iD )
!       end do !-- iD
!     end if !-- BoundaryConditions

    end select !-- PS
    end associate !-- I 

  end subroutine InitializePositionSpace


  subroutine InitializeGravitation ( FB, GravitationType )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in ) :: &
      GravitationType
!     real ( KDR ), intent ( in ), optional :: &
!       UniformAccelerationOption

    associate ( I  =>  FB % Integrator )

    select case ( trim ( GravitationType ) )
    case ( 'GALILEO' )
      allocate ( Gravitation_G_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_G_Form )
      call G % Initialize &
             ( I % X, &
               DeviceMemoryOption = FB % DeviceMemory, &
               PinnedMemoryOption = FB % PinnedMemory, &
               DevicesCommunicateOption = FB % DevicesCommunicate )
      end select !-- G
    case default
      call Show ( 'GravitationType not recognized', CONSOLE % ERROR )
      call Show ( GravitationType, 'GravitationType', CONSOLE % ERROR )
      call Show ( 'FluidBox_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GravitationType

    end associate !-- I

  end subroutine InitializeGravitation


  subroutine InitializeFluid ( FB, FluidType )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      FluidType

    select type ( I  =>  FB % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

    select case ( trim ( FluidType ) )
    case ( 'DUST' )
      allocate ( Fluid_D_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_D_Form )
      call F % Initialize ( G, FB % Units_F )
      end select !-- G
    case default
      call Show ( 'FluidType not recognized', CONSOLE % ERROR )
      call Show ( FluidType, 'FluidType', CONSOLE % ERROR )
      call Show ( 'FluidBox_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeFluid', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

    end associate !-- G
    end select !-- I

  end subroutine InitializeFluid


  subroutine InitializeStep ( FB )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB

    select type ( I  =>  FB % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )

    allocate ( Step_RK_CS_Form :: I % Step_X )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_Form )

    call S % Initialize ( F )
!     if ( present ( GravitySolverTypeOption ) ) &   
!       S % ComputeConstraints % Pointer => ComputeGravity
!       S % ApplySources % Pointer => ApplyGravity_F

    end select !-- S

    end associate !-- F
    end select !-- I

  end subroutine InitializeStep


end module FluidBox_Form
