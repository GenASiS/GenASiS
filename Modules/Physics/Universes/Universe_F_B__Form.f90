module Universe_F_B__Form

  !-- Universe_Fluid_Box__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: Universe_F_B_Form
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
  contains
    procedure, private, pass :: &
      Initialize_F_B
    generic, public :: &
      Initialize => Initialize_F_B
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_F_B
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_F_B
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeGravitation
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      InitializeStep
   end type Universe_F_B_Form


contains


  subroutine Initialize_F_B &
               ( U, FluidType, GravitationType, NameOption, &
                 MinCoordinateOption, MaxCoordinateOption, FinishTimeOption, &
                 nCellsOption, nWriteOption )
!                 CourantFactorOption, UniformAccelerationOption, &

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
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

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_F_B'

    Name  =  'Universe_F_B'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call U % Universe_H_Form % Initialize ( NameOption = Name )

    allocate ( U % Units_F ( 1 ) )

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )
    call U % InitializeGravitation &
           ( GravitationType )
!              UniformAccelerationOption = UniformAccelerationOption, &
    call U % InitializeFluid &
           ( FluidType )
    call U % InitializeStep &
           ( )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    if ( .not. allocated ( I % dT_Label ) ) then
      allocate ( I % dT_Label ( 1 ) )
      I % dT_Label ( 1 ) = 'Fluid advection'
    end if

    call I % Initialize &
           ( Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
!             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

    end select !-- I

  end subroutine Initialize_F_B


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_B_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_F ) ) &
      deallocate ( U % Units_F )

  end subroutine Finalize


  subroutine AllocateIntegrator_F_B ( U )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U

    allocate ( Integrator_CS_Form :: U % Integrator )

  end subroutine AllocateIntegrator_F_B


  subroutine InitializePositionSpace &
               ( U, MinCoordinateOption, MaxCoordinateOption, nCellsOption )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_Form  ::  I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_Form )

    call PS % Initialize &
           ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
             NameOption = 'PositionSpace', &
             CoordinateUnitOption = U % Units_F ( 1 ) % Coordinate_PS, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )

    end select !-- PS
    end associate !-- I 

  end subroutine InitializePositionSpace


  subroutine InitializeGravitation ( U, GravitationType )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
!     real ( KDR ), intent ( in ), optional :: &
!       UniformAccelerationOption

    associate ( I  =>  U % Integrator )

    select case ( trim ( GravitationType ) )
    case ( 'GALILEO' )
      allocate ( Gravitation_G_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_G_Form )
      call G % Initialize &
             ( I % X, &
               DeviceMemoryOption = U % DeviceMemory, &
               PinnedMemoryOption = U % PinnedMemory, &
               DevicesCommunicateOption = U % DevicesCommunicate )
      end select !-- G
    case default
      call Show ( 'GravitationType not recognized', CONSOLE % ERROR )
      call Show ( GravitationType, 'GravitationType', CONSOLE % ERROR )
      call Show ( 'Universe_F_B__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GravitationType

    end associate !-- I

  end subroutine InitializeGravitation


  subroutine InitializeFluid ( U, FluidType )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

    select case ( trim ( FluidType ) )
    case ( 'DUST' )
      allocate ( Fluid_D_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_D_Form )
      call F % Initialize ( G, U % Units_F )
      end select !-- G
    case ( 'IDEAL' )
      allocate ( Fluid_P_I_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_I_Form )
      call F % Initialize ( G, U % Units_F )
      end select !-- G
    case default
      call Show ( 'FluidType not recognized', CONSOLE % ERROR )
      call Show ( FluidType, 'FluidType', CONSOLE % ERROR )
      call Show ( 'Universe_F_B__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeFluid', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

    end associate !-- G
    end select !-- I

  end subroutine InitializeFluid


  subroutine InitializeStep ( U )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
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


end module Universe_F_B__Form
