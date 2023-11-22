module Universe_F_SC__Form

  !-- Universe_Fluid_SymmetricCurvilinear__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: Universe_F_SC_Form
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
  contains
    procedure, private, pass :: &
      Initialize_F_SC
    generic, public :: &
      Initialize => Initialize_F_SC
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_F_SC
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_F_SC
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeGravitation
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, public, pass :: &
      InitializeStep
  end type Universe_F_SC_Form


contains


  subroutine Initialize_F_SC &
               ( U, FluidType, Name, RadiusMax, &
                 FinishTimeOption, nCellsRadiusOption, nWriteOption )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      Name
    real ( KDR ), intent ( in ) :: &
      RadiusMax
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption, &
      nWriteOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_F_SC'

    call U % Universe_H_Form % Initialize ( Name )

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( RadiusMax, &
           nCellsRadiusOption = nCellsRadiusOption )
    call U % InitializeGravitation &
           ( )
    call U % InitializeFluid &
           ( FluidType )
    call U % SetBoundaryConditions &
           ( )
    call U % InitializeStep &
           ( )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    if ( .not. allocated ( I % dT_Label ) ) then
      allocate ( I % dT_Label ( 1 ) )
      I % dT_Label ( 1 ) = 'FluidAdvection'
    end if

    call I % Initialize &
           ( Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    end select !-- I

  end subroutine Initialize_F_SC


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_SC_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_F ) ) &
      deallocate ( U % Units_F )

  end subroutine Finalize


  subroutine AllocateIntegrator_F_SC ( U )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
      U

    allocate ( Integrator_CS_Form :: U % Integrator )

  end subroutine AllocateIntegrator_F_SC


  subroutine InitializePositionSpace ( U, RadiusMax, nCellsRadiusOption )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
      U
    real ( KDR ), intent ( in ) :: &
      RadiusMax
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_SC_Form  ::  I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_SC_Form )

    call PS % Initialize &
           ( RadiusMax, &
             CommunicatorOption = PROGRAM_HEADER % Communicator, &
             NameOption = 'PositionSpace', &
             DeviceMemoryOption = U % DeviceMemory, &
             nCellsRadiusOption = nCellsRadiusOption )

    end select !-- PS
    end associate !-- I 

  end subroutine InitializePositionSpace


  subroutine InitializeGravitation ( U )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
      U

    associate ( I  =>  U % Integrator )

    allocate ( Gravitation_G_Form  ::  I % Geometry_X )
    select type ( G  =>  I % Geometry_X )
      class is ( Gravitation_G_Form )
    call G % Initialize &
           ( I % X, &
             DeviceMemoryOption = U % DeviceMemory, &
             PinnedMemoryOption = U % PinnedMemory, &
             DevicesCommunicateOption = U % DevicesCommunicate )
    end select !-- G

    end associate !-- I

  end subroutine InitializeGravitation


  subroutine InitializeFluid ( U, FluidType )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType

    integer ( KDI ) :: &
      iB  !-- iBoundary

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

    select type ( A  =>  I % X )
    class is ( Atlas_SCG_Form )
    associate ( C  =>  A % Chart_GS )
      allocate ( U % Units_F ( 1 ) )
      call U % Units_F ( 1 ) % Initialize ( C % CoordinateUnit )
    end associate !-- C
    end select !-- A

    select case ( trim ( FluidType ) )
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
      call Show ( 'Universe_F_SC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeFluid', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

    end associate !-- G
    end select !-- I

  end subroutine InitializeFluid


  subroutine SetBoundaryConditions ( U )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    associate &
      ( F  =>  I % CurrentSet_X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_SC_Form )

    select case ( PS % Chart_GS_SC % nDimensions )
    case ( 1 )  !-- spherical coordinates

      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )

    case ( 2 )  !-- cylindrical coordinates

      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
      call F % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iC = 1, iD = 2 )

    case ( 3 )  !-- rectangular coordinates

      call F % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iC = 1, iD = 1 )
      call F % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iC = 1, iD = 2 )
      call F % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW', 'OUTFLOW' ], iC = 1, iD = 3 )

    end select !-- nDimensions

    end select !-- PS
    end associate !-- F

    end select !-- I

  end subroutine SetBoundaryConditions


  subroutine InitializeStep ( U )

    class ( Universe_F_SC_Form ), intent ( inout ) :: &
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

    DivergenceParts  =  .false.
    call PROGRAM_HEADER % GetParameter ( DivergenceParts, 'DivergenceParts' )

    if ( DivergenceParts ) then
      select type ( F )
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
        call Show ( 'Universe_F_SC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- F
    else  !-- DivergenceTotal
      select type ( F )
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
        call Show ( 'Universe_F_SC__Form', 'module', CONSOLE % ERROR )
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


end module Universe_F_SC__Form
