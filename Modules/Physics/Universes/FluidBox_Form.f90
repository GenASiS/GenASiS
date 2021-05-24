module FluidBox_Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: FluidBoxForm
    type ( Units_F_Form ), allocatable :: &
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
!     procedure, public, pass :: &
!       InitializeFluid
!     procedure, public, pass :: &
!       InitializeStep
   end type FluidBoxForm

contains


  subroutine Initialize_FB &
               ( FB, FluidType, GravitationType, Name, MinCoordinateOption, &
                 MaxCoordinateOption, nCellsOption )!, &
!                 GravitySolverTypeOption, MinCoordinateOption, &
!                 MaxCoordinateOption, FinishTimeOption, &
!                 CourantFactorOption, UniformAccelerationOption, &
!                 nWriteOption )

    class ( FluidBoxForm ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType, &
      Name
!     character ( * ), intent ( in ), optional :: &
!       GravitySolverTypeOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
!     real ( KDR ), intent ( in ), optional :: &
!       FinishTimeOption, &
!       CourantFactorOption, &
!       UniformAccelerationOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption
!     integer ( KDI ), intent ( in ), optional :: &
!       nWriteOption

    if ( FB % Type == '' ) &
      FB % Type = 'a FluidBox'

    call FB % Universe_H_Form % Initialize ( NameOption = Name )

    allocate ( FB % Units_F )

    call FB % AllocateIntegrator &
           ( )
    call FB % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )
    call FB % InitializeGravitation &
           ( GravitationType )
!              UniformAccelerationOption = UniformAccelerationOption, &
!     call FB % InitializeFluid &
!            ( FluidType, FluidUseDeviceOption = FB % UseDevice )
!     call FB % InitializeStep &
!            ( Name, GravitySolverTypeOption = GravitySolverTypeOption )


!     select type ( I => FB % Integrator )
!     class is ( Integrator_C_PS_Form )
!       call I % Initialize &
!              ( FB, Name, TimeUnitOption = FB % Units % Time, &
!                FinishTimeOption = FinishTimeOption, &
!                CourantFactorOption = CourantFactorOption, &
!                nWriteOption = nWriteOption )
!     end select !-- I

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

    allocate ( Integrator_CSA_Form :: FB % Integrator )

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
    logical ( KDL ), dimension ( 3 ) :: &
      Periodic
    
    associate ( I  =>  FB % Integrator )

    Periodic  =  .true.

    allocate ( Atlas_SCG_Form  ::  I % X_A )
    select type ( PS  =>  I % X_A )
      class is ( Atlas_SCG_Form )
    call PS % Initialize &
           ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
             NameOption = 'PositionSpace', &
             PeriodicOption = Periodic, &
             CoordinateUnitOption = FB % Units_F % Coordinate_PS, &
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
      allocate ( Gravitation_G_A_Form  ::  I % Geometry_X_A )
      select type ( GA  =>  I % Geometry_X_A )
        class is ( Gravitation_G_A_Form )
      call GA % Initialize &
             ( I % X_A, &
               DeviceMemoryOption = FB % DeviceMemory, &
               PinnedMemoryOption = FB % PinnedMemory, &
               DevicesCommunicateOption = FB % DevicesCommunicate )
      end select !-- GA
    case default
      call Show ( 'GravitationType not recognized', CONSOLE % ERROR )
      call Show ( GravitationType, 'GravitationType', CONSOLE % ERROR )
      call Show ( 'FluidBox_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    end associate !-- I

  end subroutine InitializeGravitation


!   subroutine InitializeFluid ( FB, FluidType, FluidUseDeviceOption )

!     class ( FluidBoxForm ), intent ( inout ) :: &
!       FB
!     character ( * ), intent ( in )  :: &
!       FluidType
!     logical ( KDL ), intent ( in ), optional :: &
!       FluidUseDeviceOption

!     select type ( I => FB % Integrator )
!     class is ( Integrator_C_PS_Form )

!     select type ( PS => I % PositionSpace )
!     class is ( Atlas_SC_Form )

!     allocate ( Fluid_ASC_Form :: I % Current_ASC )
!     select type ( FA => I % Current_ASC )
!     class is ( Fluid_ASC_Form )
!       call FA % Initialize &
!              ( PS, FluidType, FB % Units, &
!                UsePinnedMemoryOption = FluidUseDeviceOption )
!       if ( present ( FluidUseDeviceOption ) ) then
!         if ( FluidUseDeviceOption ) then
!           call FA % AllocateDevice ( )
!         end if
!       end if
!     end select !-- FA

!     end select !-- PS
!     end select !-- I

!   end subroutine InitializeFluid


!   subroutine InitializeStep ( FB, Name, GravitySolverTypeOption )

!     class ( FluidBoxForm ), intent ( inout ) :: &
!       FB
!     character ( * ), intent ( in )  :: &
!       Name
!     character ( * ), intent ( in ), optional :: &
!       GravitySolverTypeOption

!     select type ( I => FB % Integrator )
!     class is ( Integrator_C_PS_Form )

!     allocate ( Step_RK2_C_ASC_Form :: I % Step )
!     select type ( S => I % Step )
!     class is ( Step_RK2_C_ASC_Form )
!     call S % Initialize ( I, I % Current_ASC, NameSuffix = 'Fluid' )
!     if ( present ( GravitySolverTypeOption ) ) &   
!       S % ComputeConstraints % Pointer => ComputeGravity
!       S % ApplySources % Pointer => ApplyGravity_F
!     end select !-- S

!     end select !-- I

!   end subroutine InitializeStep


end module FluidBox_Form
