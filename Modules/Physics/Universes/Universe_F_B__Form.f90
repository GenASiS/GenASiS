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
      AllocateIntegrator
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeGravitation
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      InitializeStep
    procedure, public, pass :: &
      InitializeIntegrator
   end type Universe_F_B_Form

    private :: &
      SetSlope_N

contains


  subroutine Initialize_F_B &
               ( U, FluidType, GravitationType, Name, &
                 MinCoordinateOption, MaxCoordinateOption, FinishTimeOption, &
                 UniformAccelerationOption, nCellsOption, nWriteOption )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType, &
      Name
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      UniformAccelerationOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_F_B'

    call U % Universe_H_Form % Initialize ( Name )

    allocate ( U % Units_F ( 1 ) )
    call U % Units_F ( 1 ) % Initialize ( )

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )
    call U % InitializeGravitation &
           ( GravitationType, &
             UniformAccelerationOption = UniformAccelerationOption )
    call U % InitializeFluid &
           ( FluidType )
    call U % InitializeStep &
           ( )
    call U % InitializeIntegrator &
           ( FinishTimeOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

  end subroutine Initialize_F_B


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_B_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_F ) ) &
      deallocate ( U % Units_F )

  end subroutine Finalize


  subroutine AllocateIntegrator ( U )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U

    allocate ( Integrator_CS_Form :: U % Integrator )

  end subroutine AllocateIntegrator


  subroutine InitializePositionSpace &
               ( U, CommunicatorOption, MinCoordinateOption, &
                 MaxCoordinateOption, nCellsOption )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption

    type ( CommunicatorForm ), pointer :: &
      Communicator

    if ( present ( CommunicatorOption ) ) then
      Communicator  =>  CommunicatorOption
    else
      Communicator  =>  U % Communicator
    end if

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_Form  ::  I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_Form )

    call PS % Initialize &
           ( CommunicatorOption = Communicator, &
             NameOption = 'PositionSpace', &
             DeviceMemoryOption = U % DeviceMemory, &
             CoordinateUnitOption = U % Units_F ( 1 ) % Coordinate_PS, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )

    end select !-- PS
    end associate !-- I 

  end subroutine InitializePositionSpace


  subroutine InitializeGravitation &
               ( U, GravitationType, UniformAccelerationOption )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
    real ( KDR ), intent ( in ), optional :: &
      UniformAccelerationOption

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
    case ( 'NEWTON_UA' )

      if ( .not. present ( UniformAccelerationOption ) ) then
        call Show ( 'UniformAccelerationOption not present', CONSOLE % ERROR )
        call Show ( 'NEWTON_UA', 'GravitationType', CONSOLE % ERROR )
        call Show ( 'Universe_F_B__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      allocate ( Gravitation_N_UA_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_N_UA_Form )
      call G % Initialize &
             ( I % X, &
               Acceleration = UniformAccelerationOption, &
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

    integer ( KDI ) :: &
      iB  !-- iBoundary

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

    select case ( trim ( FluidType ) )
    case ( 'DUST' )

      allocate ( Fluid_D_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )

        !-- TallyInterior must be allocated before F % Initialize...
        allocate ( Tally_F_D_Form :: F % TallyInterior )
        allocate ( Tally_F_D_Form :: F % TallyTotal )
        allocate ( Tally_F_D_Form :: F % TallyChange )
        select type ( TI  =>  F % TallyInterior )
          type is ( Tally_F_D_Form )
        call TI % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TI
        select type ( TT  =>  F % TallyTotal )
          type is ( Tally_F_D_Form )
        call TT % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TT
        select type ( TC  =>  F % TallyChange )
          type is ( Tally_F_D_Form )
        call TC % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TC

        !-- ( Initialize the Fluid)
        call F % Initialize ( G, U % Units_F )

        !-- ... but TallyBoundary needs F % Initialize already called.
        allocate ( F % TallyBoundary ( F % nBoundaries ) )
        do iB  =  1,  F % nBoundaries
          allocate ( Tally_F_D_Form :: F % TallyBoundary ( iB ) % Element )
          select type ( TB  =>  F % TallyBoundary ( iB ) % Element )
            type is ( Tally_F_D_Form )
          call TB % Initialize ( G, U % Units_F ( 1 ) )
          end select !-- TB
        end do !-- iB

        !-- Boundary accumulation storage
        call F % AllocateBoundary_SCG ( nT = F % TallyInterior % nSelected )

      end select !-- F

    case ( 'IDEAL' )

      allocate ( Fluid_P_I_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_I_Form )

        !-- TallyInterior must be allocated before F % Initialize...
        allocate ( Tally_F_P_Form :: F % TallyInterior )
        allocate ( Tally_F_P_Form :: F % TallyTotal )
        allocate ( Tally_F_P_Form :: F % TallyChange )
        select type ( TI  =>  F % TallyInterior )
          type is ( Tally_F_P_Form )
        call TI % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TI
        select type ( TT  =>  F % TallyTotal )
          type is ( Tally_F_P_Form )
        call TT % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TT
        select type ( TC  =>  F % TallyChange )
          type is ( Tally_F_P_Form )
        call TC % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TC

        !-- ( Initialize the Fluid )
        call F % Initialize ( G, U % Units_F )

        !-- ... but TallyBoundary needs F % Initialize already called.
        allocate ( F % TallyBoundary ( F % nBoundaries ) )
        do iB  =  1,  F % nBoundaries
          allocate ( Tally_F_P_Form :: F % TallyBoundary ( iB ) % Element )
          select type ( TB  =>  F % TallyBoundary ( iB ) % Element )
            type is ( Tally_F_P_Form )
          call TB % Initialize ( G, U % Units_F ( 1 ) )
          end select !-- TB
        end do !-- iB

        !-- Boundary accumulation storage
        call F % AllocateBoundary_SCG ( nT = F % TallyInterior % nSelected )

      end select !-- F

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

    logical ( KDL ) :: &
      DivergenceParts
    character ( LDL ) :: &
      RiemannSolverType

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    associate &
      ( G  =>  I % Geometry_X )

    allocate ( Step_RK_CS_Form :: I % Step_X )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_Form )

    select type ( G )
    class is ( Gravitation_N_H_Form )
      S % SetSlope  =>  SetSlope_N
    end select !-- G

    DivergenceParts  =  .false.
    call PROGRAM_HEADER % GetParameter ( DivergenceParts, 'DivergenceParts' )

    if ( DivergenceParts ) then
      select type ( F )
      class is ( Fluid_D_Form )
        allocate ( S % DivergencePart ( 1 ) )
        associate ( DP_1D  =>  S % DivergencePart )
          allocate ( DivergencePart_F_D_V_Form :: DP_1D ( 1 ) % Element )
          associate ( DP  =>  DP_1D ( 1 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
        end associate !-- DP_1D
      class is ( Fluid_P_I_Form )
        allocate ( S % DivergencePart ( 2 ) )
        associate ( DP_1D  =>  S % DivergencePart )
          allocate ( DivergencePart_F_P_V_Form :: DP_1D ( 1 ) % Element )
          associate ( DP  =>  DP_1D ( 1 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
          allocate ( DivergencePart_F_P_P_Form :: DP_1D ( 2 ) % Element )
          associate ( DP  =>  DP_1D ( 2 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
        end associate !-- DP_1D
      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_F_B__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- F
    else  !-- DivergenceTotal
      select type ( F )
      class is ( Fluid_D_Form )
        allocate ( DivergencePart_F_D_T_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
          call DT % Initialize ( F )
        end associate !-- DT
      class is ( Fluid_P_I_Form )

        allocate ( DivergencePart_F_P_T_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
          call DT % Initialize ( F )
        end associate !-- DT

        RiemannSolverType = 'HLLC'
        call PROGRAM_HEADER % GetParameter &
               ( RiemannSolverType, 'RiemannSolverType' )
        if ( trim ( RiemannSolverType ) == 'HLLC' ) then
          allocate ( RiemannSolver_HLLC_P_Form :: S % RiemannSolver )
          associate ( RS  =>  S % RiemannSolver )
          call RS % Initialize ( F )
          end associate !-- RS
        end if
        
      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_F_B__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- F
    end if  !-- DivergenceParts

    call S % Initialize ( F )

    end select !-- S

    end associate !-- G
    end associate !-- F
    end select !-- I

  end subroutine InitializeStep


  subroutine InitializeIntegrator ( U, FinishTimeOption, nWriteOption )

    class ( Universe_F_B_Form ), intent ( inout ) :: &
      U
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    if ( .not. allocated ( I % dT_Label ) ) then
      allocate ( I % dT_Label ( 1 ) )
      I % dT_Label ( 1 ) = 'FluidAdvection'
    end if

    call I % Initialize &
           ( CommunicatorOption = U % Communicator, &
             Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    end select !-- I

  end subroutine InitializeIntegrator


  subroutine SetSlope_N ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    character ( 1 ) :: &
      StageNumber

    allocate ( Slope_DFV_N_Form :: K )

    select type ( K )
      class is ( Slope_DFV_N_Form )
    select type ( S )
      class is ( Step_RK_CS_Form )
    select type ( F  =>  S % CurrentSet )
      class is ( Fluid_D_Form ) 

    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_B ( 1 ) )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_B ( 2 ) )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_B ( 3 ) )

    !-- Dust
    iEnergy_B  =  0

    !-- Perfect fluid
    select type ( F )
    class is ( Fluid_P_Form )
      call Search &
             ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_B )
    end select !-- F

    if ( allocated ( S % DivergenceTotal ) ) then
      call K % Initialize &
             ( S % RiemannSolver, &
               S % DivergenceTotal, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B )
    else if ( allocated ( S % DivergencePart ) ) then
      call K % Initialize &
             ( S % RiemannSolver, &
               S % DivergencePart, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B )
    end if  !-- DivergenceTotal

    end select !-- F
    end select !-- S
    end select !-- K

  end subroutine SetSlope_N


end module Universe_F_B__Form
