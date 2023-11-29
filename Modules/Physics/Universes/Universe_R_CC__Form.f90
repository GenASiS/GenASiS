module Universe_R_CC__Form

  !-- Universe_Radiation_CentralCore__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Universe_F_CC__Form

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: Universe_R_CC_Form
    integer ( KDI ) :: &
      iRadiation  = 0, &
      nRadiations = 0
    real ( KDR ) :: &
      InteractionFactor = 0.0_KDR
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
    type ( CommunicatorForm ), allocatable :: &
      Communicator_PS  !-- PositionSpace
!    type ( CollectiveOperation_R_Form ), dimension ( : ), allocatable :: &
!      CO_SplitSource
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
    class ( Interactions_NM_G_Form ), allocatable :: &
      Interactions_NM_G
  contains
    procedure, private, pass :: &
      Initialize_R_CC
    generic, public :: &
      Initialize => Initialize_R_CC
    final :: &
      Finalize
    procedure, private, pass :: &
      SetCommunicator
    procedure, private, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      InitializeInteractions
    procedure, public, pass :: &
      SetBoundaryConditions
!    procedure, public, pass :: &
!      InitializeSteps
    procedure, public, pass :: &
      InitializeStep
    procedure, public, pass :: &
      InitializeIntegrator
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass ( U ) :: &
      Compute_dT_RT_CGS
  end type Universe_R_CC_Form

    !-- FIXME: This is for a workaround in SetSlope routines below
    class ( Universe_R_CC_Form ), private, pointer :: &
      UNIVERSE => null ( )

    private :: &
!      ResolveCycle_R, &
!      PrepareStep_F, &
!      ComputeSource_F, &
      Compute_dT_Local, &
      InitializeSeries, &
      Analyze, &
      Set_T_CheckpointInterval, &
      SetSlope_F_P_DFV_N_SS, &
      SetSlope_NM_G_I, &
      SetSlope_NM_G_DFV_I

      private :: &
        Compute_dT_RT_CGS_Kernel

    interface
    
      module subroutine Compute_dT_RT_CGS_Kernel &
               ( dT_E, dT_N, ProperCell, Q, R, E, N, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT_E, dT_N
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Q, R, &
          E, N
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_RT_CGS_Kernel

    end interface


contains


  subroutine Initialize_R_CC &
               ( U, RadiationName, RadiationType, FormalismType, FluidType, &
                 GravitationType, Name, UnitsTypeOption, FinishTimeOption, &
                 RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_R_CC_Form ), intent ( inout ), target :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      FluidType, &
      GravitationType, &
      Name
    character ( * ), intent ( in ), optional :: &
      UnitsTypeOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_CC'

    call U % Universe_H_Form % Initialize &
           ( Name, UnitsTypeOption = UnitsTypeOption )

    !-- FIXME: This is for a workaround in SetSlope routines below
    UNIVERSE  =>  U

    !-- Radiations

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    if ( trim ( RadiationType ( 1 ) )  ==  'NONE' ) then
      call U % Universe_F_CC_Form % Initialize &
             ( FluidType, GravitationType, Name, &
               FinishTimeOption = FinishTimeOption, &
               RadiusMaxOption = RadiusMaxOption, &
               RadiusCoreOption = RadiusCoreOption, &
               RadialRatioOption = RadialRatioOption, &
               nCellsPolarOption = nCellsPolarOption, &
               nWriteOption = nWriteOption )
      return
    end if

    !-- Initializations

    call U % SetCommunicator &
           ( )
    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( CommunicatorOption = U % Communicator_PS, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )
    ! call RB % InitializeMomentumSpace &
    !        ( EnergySpacingOption = EnergySpacingOption, &
    !          MinEnergyOption = MinEnergyOption, &
    !          MaxEnergyOption = MaxEnergyOption, &
    !          MinWidthEnergyOption = MinWidthEnergyOption, &
    !          EnergyScaleOption = EnergyScaleOption, &
    !          nCellsEnergyOption = nCellsEnergyOption )
    call U % InitializeGravitation &
           ( GravitationType )
    call U % InitializeFluid &
           ( FluidType )
    call U % InitializeRadiation &
           ( )
    call U % InitializeInteractions &
           ( )
    call U % SetBoundaryConditions &
           ( )
!    call U % InitializeSteps &
!           ( )
    call U % InitializeStep &
           ( )
    call U % InitializeIntegrator &
           ( GravitationType, &
             FinishTimeOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    call U % SetMeasures ( )

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
!    I % ResolveCycle              =>  ResolveCycle_R
!    I % PrepareStep               =>  PrepareStep_F
    I % Compute_dT_Local          =>  Compute_dT_Local
    I % InitializeSeries          =>  InitializeSeries
    I % Analyze                   =>  Analyze
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    end associate !-- I

  end subroutine Initialize_R_CC


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Interactions_NM_G ) ) &
      deallocate ( U % Interactions_NM_G )
    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
!    if ( allocated ( U % CO_SplitSource ) ) &
!      deallocate ( U % CO_SplitSource )
    if ( allocated ( U % Communicator_PS ) ) &
      deallocate ( U % Communicator_PS )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine SetCommunicator ( U )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iPS, &  !-- iPositionSpace
      nProcesses, &
      nProcesses_PS  !-- PositionSpace
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank

    nProcesses  =  U % Communicator % Size
    if ( mod ( nProcesses, U % nRadiations )  /=  0 ) then
      call Show ( 'nRadiations must evenly divide nProcesses', CONSOLE % ERROR )
      call Show ( nProcesses, 'nProcesses', CONSOLE % ERROR )
      call Show ( U % nRadiations, 'nRadiations', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )

      nProcesses_PS  =  nProcesses / U % nRadiations

      do iPS  =  1,  U % nRadiations
        allocate ( Rank, &
                   source =  [ ( iP, iP = ( iPS - 1 ) * nProcesses_PS, &
                                            iPS * nProcesses_PS  -  1 ) ] )
        if ( any ( U % Communicator % Rank  ==  Rank ) ) then
          U % iRadiation = iPS
          allocate ( U % Communicator_PS )
          call U % Communicator_PS % Initialize &
                 ( U % Communicator, Rank, NameOption = 'Communicator_PS' )
        end if
        deallocate ( Rank )
      end do !-- iPS

    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCommunicator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine SetCommunicator


  subroutine AllocateIntegrator ( U )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )
      allocate ( Integrator_CS_1D_BM_CS_Form :: U % Integrator )
      select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
        allocate ( I % Communicator_X_1D )
      end select !-- I
    case ( 'SPECTRAL' )
      allocate ( Integrator_CS_1D_CB_CS_Form :: U % Integrator )
    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine AllocateIntegrator


  subroutine InitializeRadiation ( U )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )

      select type ( A  =>  I % X )
      class is ( Atlas_SCG_Form )
      associate ( C  =>  A % Chart_GS )
        allocate ( U % Units_R ( 1 ) )
        call U % Units_R ( 1 ) % Initialize &
               ( C % CoordinateUnit, TypeOption = U % UnitsType )
      end associate !-- C
      end select !-- A

      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_Form )
      associate &
        ( iR  =>  U % iRadiation )

      select case ( trim ( U % RadiationType ( iR ) ) )
      case ( 'NEUTRINOS_E', 'NEUTRINOS_E_BAR' )

        allocate ( NeutrinoMoments_G_Form :: I % CurrentSet_X_1D )
        select type ( R  =>  I % CurrentSet_X_1D )
        class is ( NeutrinoMoments_G_Form )

        call R % Initialize &
               ( F, U % Units_R, U % RadiationType ( iR ), &
                 NameOption = U % RadiationName ( iR ) )
        if ( allocated ( U % Interactions_NM_G ) ) &
          call R % SetInteractions ( U % Interactions_NM_G )

!     select type ( I  =>  R % Interactions )
!     class is ( Interactions_BM_Form )
! call Show ( '>>> Interactions_BM InitializeRadiation' )
!     end select

!     select type ( I  =>  R % Interactions )
!     type is ( Interactions_NM_G_Form )
! call Show ( '>>> Interactions_NM_G InitializeRadiation' )
!     end select

        end select !-- R

      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( U % RadiationType ( iR ), 'RadiationType', &
                    CONSOLE % ERROR )
        call Show ( 'Universe_R_CC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- RadiationType

      end associate !-- iR
      end select !-- F

    end select !-- I

  end subroutine InitializeRadiation


  subroutine SetBoundaryConditions ( U )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U
    
    call U % Universe_F_CC_Form % SetBoundaryConditions ( )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )

    associate &
      ( F  =>  I % CurrentSet_X_1D )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    call F % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end associate !-- F

    end select !-- I

  end subroutine SetBoundaryConditions


  subroutine InitializeInteractions ( U )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( RadiationMoments_BM_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )

    if ( allocated ( U % Interactions_NM_G ) ) &
      call U % Interactions_NM_G % Initialize ( R, U % Units_R, F )

    end select !-- F
    end select !-- R
    end select !-- I

  end subroutine InitializeInteractions


!   subroutine InitializeSteps ( U )

!     class ( Universe_R_CC_Form ), intent ( inout ) :: &
!       U

!     integer ( KDI ) :: &
!       EvolutionOrder
!     character ( LDL ) :: &
!       RiemannSolverType

!     EvolutionOrder  =  2
!     call PROGRAM_HEADER % GetParameter ( EvolutionOrder, 'EvolutionOrder' )

!     U % Coarsen  =  .true.
!     call PROGRAM_HEADER % GetParameter ( U % Coarsen, 'Coarsen' )

!     !-- Fluid step

!     select type ( I  =>  U % Integrator )
!       class is ( Integrator_CS_Form )
!     select type ( F  =>  I % CurrentSet_X )
!       class is ( Fluid_P_HN_Form )

!     allocate ( Step_RK_CS_Form :: I % Step_X )
!     select type ( S  =>  I % Step_X )
!       class is ( Step_RK_CS_Form )

!     allocate ( DivergencePart_F_P_HN_T_Form :: S % DivergenceTotal )
!     associate ( DT  =>  S % DivergenceTotal )
!       call DT % Initialize ( F )
!     end associate !-- DT

!     RiemannSolverType = 'HLL'
!     call PROGRAM_HEADER % GetParameter &
!            ( RiemannSolverType, 'RiemannSolverType' )
!     if ( trim ( RiemannSolverType ) == 'HLLC' ) then
!       allocate ( RiemannSolver_HLLC_P_HN_Form :: S % RiemannSolver )
!       associate ( RS  =>  S % RiemannSolver )
!       call RS % Initialize ( F )
!       end associate !-- RS
!     end if        

!     S % SetSlope  =>  SetSlope_F_P_DFV_N_SS

!     call S % Initialize ( F, OrderOption = EvolutionOrder )

!     !-- Coarsening
!     if ( U % Coarsen ) then
!       allocate ( U % Coarsening )
!       associate &
!         ( C  =>  U % Coarsening, &
!           G  =>  I % Geometry_X )
!         call C % Initialize ( F, G )
!         call S % SetCoarsening ( C )
!       end associate !-- C, etc.
!     end if

!     end select !-- S

!     class default
!       call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
!       call Show ( 'Universe_R_CC__Form', 'module', CONSOLE % ERROR )
!       call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
!       call PROGRAM_HEADER % Abort ( )
!     end select    !-- F
!     end select !-- I

!     !-- Radiation step

!     select type ( I  =>  U % Integrator )
!     class is ( Integrator_CS_1D_BM_CS_Form )

!       associate &
!         ( R  =>  I % CurrentSet_X_1D )

!       allocate ( Step_RK_CS_Form :: I % Step_1D )
!       select type ( S  =>  I % Step_1D )
!         class is ( Step_RK_CS_Form )

!       allocate ( DivergencePart_NM_G_Form :: S % DivergenceTotal )
!       associate ( DT  =>  S % DivergenceTotal )
!       call DT % Initialize ( R )
!       end associate !-- DT

!       allocate ( DiffusionFactor_RM_Form :: S % DiffusionFactor )
!       select type ( DF  =>  S % DiffusionFactor )
!       class is ( DiffusionFactor_RM_Form )
!         call DF % Initialize ( U % Interactions_NM_G )
!       end select !-- DF

! !      S % SetSlope  =>  SetSlope_NM_G_I
!       S % SetSlope  =>  SetSlope_NM_G_DFV_I

!       call S % Initialize ( R, OrderOption = EvolutionOrder )

!       end select !-- S
!       end associate !-- R

!     end select !-- I

!   end subroutine InitializeSteps


  subroutine InitializeStep ( U )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      EvolutionOrder
    character ( LDL ) :: &
      RiemannSolverType

    EvolutionOrder  =  2
    call PROGRAM_HEADER % GetParameter ( EvolutionOrder, 'EvolutionOrder' )

    U % Coarsen  =  .true.
    call PROGRAM_HEADER % GetParameter ( U % Coarsen, 'Coarsen' )

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )

      allocate ( Step_RK_CS_CS_Form :: I % Step_X )
      select type ( S  =>  I % Step_X )
        class is ( Step_RK_CS_CS_Form )

      allocate ( S % Step_CS_1 )
      allocate ( S % Step_CS_2 )
      associate &
        ( S_R  =>  S % Step_CS_1, &
          S_F  =>  S % Step_CS_2 )

      !-- Radiation
      associate ( R  =>  I % CurrentSet_X_1D )

      allocate ( DivergencePart_NM_G_Form :: S_R % DivergenceTotal )
      associate ( DT  =>  S_R % DivergenceTotal )
      call DT % Initialize ( R )
      end associate !-- DT

      allocate ( DiffusionFactor_RM_Form :: S_R % DiffusionFactor )
      select type ( DF  =>  S_R % DiffusionFactor )
      class is ( DiffusionFactor_RM_Form )
        call DF % Initialize ( U % Interactions_NM_G )
      end select !-- DF

      S_R % SetSlope  =>  SetSlope_NM_G_DFV_I

      !-- Fluid
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_HN_Form )

      allocate ( DivergencePart_F_P_HN_T_Form :: S_F % DivergenceTotal )
      associate ( DT  =>  S_F % DivergenceTotal )
        call DT % Initialize ( F )
      end associate !-- DT

      RiemannSolverType = 'HLL'
      call PROGRAM_HEADER % GetParameter &
             ( RiemannSolverType, 'RiemannSolverType' )
      if ( trim ( RiemannSolverType ) == 'HLLC' ) then
        allocate ( RiemannSolver_HLLC_P_HN_Form :: S_F % RiemannSolver )
        associate ( RS  =>  S_F % RiemannSolver )
        call RS % Initialize ( F )
        end associate !-- RS
      end if        

      S_F % SetSlope  =>  SetSlope_F_P_DFV_N_SS

      !-- Combined step
      call S % Initialize ( R, F, OrderOption = EvolutionOrder )

      !-- Coarsening
      if ( U % Coarsen ) then
        allocate ( U % Coarsening )
        associate &
          ( C  =>  U % Coarsening, &
            G  =>  I % Geometry_X )
          call C % Initialize ( F, G )
          call S_F % SetCoarsening ( C )
        end associate !-- C, etc.
      end if

      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_R_CC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select    !-- F
      end associate !-- R
      end associate !-- S_R, S_F
      end select !-- S

    end select !-- I

  end subroutine InitializeStep


  subroutine InitializeIntegrator &
               ( U, GravitationType, FinishTimeOption, GravityFactorOption, &
                 nWriteOption )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      GravityFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iCS  !-- iCurrentSet

    select type ( I => U % Integrator )
    class is ( Integrator_CS_1D_CS_Form )

    I % iCurrentSet  =  U % iRadiation

    I % nCurrentSets  =  size ( U % RadiationName )
    allocate ( I % dT_Label ( 5 ) )

    I % dT_Label ( 1 )  =  'GravitationAcceleration'
    I % dT_Label ( 2 )  =  'FluidAdvection'
    I % dT_Label ( 3 )  =  'RadiationStreaming'
    I % dT_Label ( 4 )  =  'EnergyTransfer'
    I % dT_Label ( 5 )  =  'ElectronNumberTransfer'

    U % GravityFactor  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter &
           ( U % GravityFactor, 'GravityFactor' )

    U % InteractionFactor  =  1.0e-2_KDR
    call PROGRAM_HEADER % GetParameter &
           ( U % InteractionFactor, 'InteractionFactor' )

    I % StreamSuffix  =  '_' // trim ( U % RadiationName ( U % iRadiation ) )

    call I % Initialize &
           ( CommunicatorOption = U % Communicator_PS, &
             Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    call U % Interactions_NM_G % SetStream ( I % Checkpoint_X )
    call U % Interactions_NM_G % Show ( )

    end select !-- I

  end subroutine InitializeIntegrator


  subroutine ShowParameters ( U )

    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U

    call U % Universe_F_CC_Form % ShowParameters ( )

    call Show ( U % RadiationName,     'RadiationName',     U % IGNORABILITY )
    call Show ( U % RadiationType,     'RadiationType',     U % IGNORABILITY )
    call Show ( U % iRadiation,        'iRadiation',        U % IGNORABILITY )
    call Show ( U % FormalismType,     'FormalismType',     U % IGNORABILITY )
    call Show ( U % InteractionFactor, 'InteractionFactor', U % IGNORABILITY )

  end subroutine ShowParameters


  subroutine Compute_dT_RT_CGS ( dT_E, dT_N, U, iC, T_Option )

    real ( KDR ), intent ( inout ) :: &
      dT_E, dT_N
    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iEnergy_B, iNumber_B

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_CS_Form )
    select type ( S_R  =>  S % Step_CS_1 )
      class is ( Step_RK_CS_Form )
    select type ( S_R_I  =>  S_R % Slope % Component ( 2 ) % Element )
      class is ( Slope_NM_G_I_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_HN_Form )
    associate &
      ( FS =>  F % SplitSource )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (   C  =>   A % Chart_GS, &
        FSV  =>  FS % Storage_GS % Value, &
         FV  =>   F % Storage_GS % Value )

    call Search ( F % iaBalanced, F % ENERGY_DENSITY_B,   iEnergy_B )
    call Search ( F % iaBalanced, F % ELECTRON_DENSITY_B, iNumber_B )

    call Compute_dT_RT_CGS_Kernel &
           ( dT_E, dT_N, C % ProperCell, &
             Q  =  FSV ( :, iEnergy_B ), &
             R  =  FSV ( :, iNumber_B ), &
             E  =   FV ( :, F % ENERGY_DENSITY_B ), &
             N  =   FV ( :, F % ELECTRON_DENSITY_B ), &
             UseDeviceOption = F % DeviceMemory )

    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Universe_R_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_dT_RT_CGS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- FS
    end select !-- F
    end select !-- S_R_I
    end select !-- S_R
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_RT_CGS


!   subroutine ResolveCycle_R ( I )

!     class ( Integrator_H_Form ), intent ( inout ) :: &
!       I

!     select type ( U  =>  I % System )
!       class is ( Universe_R_CC_Form )
!     select type ( I )
!       class is ( Integrator_CS_1D_BM_CS_Form )
!     select type ( R  =>  I % CurrentSet_X_1D )
!       class is ( NeutrinoMoments_G_Form )

! !    call R % SetFluidVelocity ( )
! !    call R % ComputeFromBalanced ( )
!     call R % ComputeSpectralParameters ( )
!     call R % ComputeEquilibrium ( )
!     call U % Interactions_NM_G % Compute ( )

!     select type ( S_1D  =>  I % Step_1D )
!       class is ( Step_RK_CS_Form )
!     if ( S_1D % Slope % nComponents  >  1 ) then
!       select type ( S_R_I  =>  S_1D % Slope % Component ( 2 ) % Element )
!         class is ( Slope_NM_G_I_Form )

!       !-- To be used for EnergyTransfer and ElectronNumberTransfer time steps
!       call S_R_I % Compute ( dT = 0.0_KDR )
!       call ComputeSource_F ( I, S_R_I )

!       end select !-- S_R_I
!     end if !-- Slope % nComponents > 1
!     end select !-- S_1D

!     end select !-- R
!     end select !-- I
!     end select !-- U

!   end subroutine ResolveCycle_R


!   subroutine PrepareStep_F ( I )

!     class ( Integrator_H_Form ), intent ( inout ) :: &
!       I

!     select type ( I )
!       class is ( Integrator_CS_1D_BM_CS_Form )
!     select type ( S_1D  =>  I % Step_1D )
!       class is ( Step_RK_CS_Form )
!     if ( S_1D % Slope % nComponents  >  1 ) then
!       select type ( SS_R_I  =>  S_1D % SlopeSum % Component ( 2 ) % Element )
!         class is ( Slope_NM_G_I_Form )

!       call ComputeSource_F ( I, SS_R_I )

!       end select !-- SS_R_I
!     end if !-- Slope % nComponents > 1
!     end select !-- S_1D
!     end select !-- I

!   end subroutine PrepareStep_F


!   subroutine ComputeSource_F ( I, S_R_I )

!     class ( Integrator_H_Form ), intent ( inout ) :: &
!       I
!     class ( Slope_NM_G_I_Form ), intent ( inout ) :: &
!       S_R_I

!     integer ( KDI ) :: &
!       iC, &  !-- iChart
!       iEnergy_R, iEnergy_F, &
!       iNumber_R, iNumber_F, &
!       nSources, &
!       nValues
!     integer ( KDI ), dimension ( 3 ) :: &
!       iMomentum_R, iMomentum_F
!     real ( KDR ) :: &
!       NumberFactor
!     real ( KDR ), dimension ( :, : ), pointer :: &
!       RSB, &  !-- 2D alias for outgoing buffer
!       FSB     !-- 2D alias for incoming buffer

!     select type ( U  =>  I % System )
!       class is ( Universe_R_CC_Form )
!     select type ( I )
!       class is ( Integrator_CS_1D_BM_CS_Form )
!     select type ( R  =>  I % CurrentSet_X_1D )
!       class is ( NeutrinoMoments_G_Form )
!     select type ( F  =>  I % CurrentSet_X )
!       class is ( Fluid_P_HN_Form )

!     call Search &
!            ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
!     call Search &
!            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
!     call Search &
!            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
!     call Search &
!            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )
!     call Search &
!            ( R % iaBalanced, R % NUMBER_DENSITY_B, iNumber_R )

!     call Search &
!            ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
!     call Search &
!            ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
!     call Search &
!            ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
!     call Search &
!            ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )
!     call Search &
!            ( F % iaBalanced, F % ELECTRON_DENSITY_B, iNumber_F )

!     nSources  =  5

!     if ( .not. allocated ( U % CO_SplitSource ) ) &
!       allocate ( U % CO_SplitSource ( F % Atlas % nCharts ) )

!     do iC  =  1,  F % Atlas % nCharts
!       associate &
!         ( CO  =>  U % CO_SplitSource ( iC ) )
!       associate &
!         ( RSV  =>  S_R_I % Storage ( iC ) % Value, &
!           FSV  =>  F % SplitSource % Storage ( iC ) % Value )
!       associate &
!         ( RS_E    =>  RSV ( :, iEnergy_R ), &
!           RS_S_1  =>  RSV ( :, iMomentum_R ( 1 ) ), &
!           RS_S_2  =>  RSV ( :, iMomentum_R ( 2 ) ), &
!           RS_S_3  =>  RSV ( :, iMomentum_R ( 3 ) ), &
!           RS_D    =>  RSV ( :, iNumber_R ), &
!           FS_G    =>  FSV ( :, iEnergy_F ), &
!           FS_S_1  =>  FSV ( :, iMomentum_F ( 1 ) ), &
!           FS_S_2  =>  FSV ( :, iMomentum_F ( 2 ) ), &
!           FS_S_3  =>  FSV ( :, iMomentum_F ( 3 ) ), &
!           FS_D    =>  FSV ( :, iNumber_F ) )

!       nValues  =  size ( FSV, dim = 1 )

!       if ( .not. allocated ( CO % Outgoing ) ) then
!         call CO % Initialize &
!                ( I % Communicator_X_1D, &
!                  nOutgoing  =  [ nValues * nSources ], &
!                  nIncoming  =  [ nValues * nSources ] )
!         if ( S_R_I % DeviceMemory .and. S_R_I % DevicesCommunicate ) then
!           call CO % AllocateDevice ( )
!         end if
!       end if
      
!       if ( .not. CO % AllocatedDevice ) &
!         call S_R_I % UpdateHost ( )
      
!       RSB ( 1 : nValues,  1 : nSources )  =>  CO % Outgoing % Value
!       FSB ( 1 : nValues,  1 : nSources )  =>  CO % Incoming % Value

!       call Copy ( RS_E,   RSB ( :, 1 ), &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( RS_S_1, RSB ( :, 2 ), &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( RS_S_2, RSB ( :, 3 ), &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( RS_S_3, RSB ( :, 4 ), &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( RS_D,   RSB ( :, 5 ), &
!                   UseDeviceOption = CO % AllocatedDevice )
      
!       !-- Energy / Momentum
!       call Multiply ( RSB ( :, 1 : 4 ), -1.0_KDR, &
!                       UseDeviceOption = CO % AllocatedDevice )

!       !-- Electron number
!       select case ( trim ( R % RadiationType ) )
!       case ( 'NEUTRINOS_E' )
!         NumberFactor  =  - 1.0_KDR
!       case ( 'NEUTRINOS_E_BAR' )
!         NumberFactor  =  + 1.0_KDR
!       case default
!         NumberFactor  =    0.0_KDR
!       end select !-- RadiationType
!       call Multiply ( RSB ( :, 5 ), NumberFactor, &
!                       UseDeviceOption = CO % AllocatedDevice )

!       call CO % Reduce ( REDUCTION % SUM )

!       call Copy ( FSB ( :, 1 ), FS_G, &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( FSB ( :, 2 ), FS_S_1, &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( FSB ( :, 3 ), FS_S_2, &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( FSB ( :, 4 ), FS_S_3, &
!                   UseDeviceOption = CO % AllocatedDevice )
!       call Copy ( FSB ( :, 5 ), FS_D, &
!                   UseDeviceOption = CO % AllocatedDevice )
      
!       if ( .not. CO % AllocatedDevice ) &
!         call F % SplitSource % UpdateDevice ( )

!       end associate !-- RS_E, etc.
!       end associate !-- RSV, etc.
!       end associate !-- CO
!     end do !-- iC

!     end select !-- F
!     end select !-- R
!     end select !-- I
!     end select !-- U

!   end subroutine ComputeSource_F


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_Form )
    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
    associate &
      ( dT_1  =>  dT_Candidate ( 1 ), &
        dT_2  =>  dT_Candidate ( 2 ), &
        dT_3  =>  dT_Candidate ( 3 ), &
        dT_4  =>  dT_Candidate ( 4 ), &
        dT_5  =>  dT_Candidate ( 5 ) )

    !-- Gravity step

    call U % Compute_dT_G_CGS ( dT_1, iC, T_Option )
    dT_1  =  U % GravityFactor  *  dT_1    

    !-- Fluid advection step

    if ( U % Coarsen ) then
      call U % Compute_dT_CS_CGS_C ( dT_2, iC, T_Option )
    else !-- .not. Coarsen
      call I % Compute_dT_CS_CGS &
             ( I % EigenspeedSet_X, dT_2, iC, T_Option )
    end if !-- Coarsen
    dT_2  =  I % CourantFactor  *  dT_2
    
    !-- Radiation streaming step

    call I % Compute_dT_CS_CGS &
           ( I % EigenspeedSet_X_1D, dT_3, iC, T_Option )
    dT_3  =  I % CourantFactor_1D  *  dT_3

    !-- Radiative transfer steps

    call U % Compute_dT_RT_CGS ( dT_4, dT_5, iC, T_Option )
    dT_4  =  U % InteractionFactor  *  dT_4
    dT_5  =  U % InteractionFactor  *  dT_5

    !-- Reduce across radiation types

    call CO % Initialize &
           ( I % Communicator_X_1D, nOutgoing = [ 3 ], &
             nIncoming = [ 3 ] )

    CO % Outgoing % Value  =  I % dT_Candidate ( 3 : 5 )
    call CO % Reduce ( REDUCTION % MIN )
    I % dT_Candidate ( 3 : 5 )  =  CO % Incoming % Value

    end associate !-- dT_1, etc.
    end select !-- I
    end select !-- U

  end subroutine Compute_dT_Local


  subroutine InitializeSeries ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_Form )

    call U % InitializeSeries_F_CC ( I )

    end select !-- U

  end subroutine InitializeSeries


  subroutine Analyze ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_Form )

    call U % Analyze_F_CC ( I, Ignorability, T_Option )

    end select !-- U

  end subroutine Analyze


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_Form )

    call U % Set_T_CheckpointInterval_F_CC ( I )

    end select !-- U

  end subroutine Set_T_CheckpointInterval


  subroutine SetSlope_F_P_DFV_N_SS ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_F_P_DFV_N_SS_Form :: K )
    select type ( K )
      class is ( Slope_F_P_DFV_N_SS_Form )
    select type ( F  =>  S % CurrentSet )
      class is ( Fluid_P_Form )

    call K % Initialize &
           ( S % RiemannSolver, S % DiffusionFactor, S % DivergenceTotal, F )
!             IgnorabilityOption = S % IGNORABILITY )

    end select !-- F
    end select !-- K
    end select !-- S

  end subroutine SetSlope_F_P_DFV_N_SS


  subroutine SetSlope_NM_G_I ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_NM_G_I_Form :: K )
    select type ( K )
      class is ( Slope_NM_G_I_Form )
    select type ( R  =>  S % CurrentSet )
      class is ( NeutrinoMoments_G_Form )

!     select type ( I  =>  R % Interactions )
!     class is ( Interactions_BM_Form )
! call Show ( '>>> Interactions_BM SetSlope_NM_G_I' )
!     end select

!     select type ( I  =>  R % Interactions )
!     type is ( Interactions_NM_G_Form )
! call Show ( '>>> Interactions_NM_G SetSlope_NM_G_I' )
!     end select

!     select type ( I  =>  UNIVERSE % Interactions_NM_G )
!     class is ( Interactions_NM_G_Form )
! call Show ( '>>> Interactions_NM_G SetSlope_NM_G_I UNIVERSE' )
!     end select

    call K % Initialize &
           ( R )!, &
!             IgnorabilityOption = S % IGNORABILITY )

    !-- FIXME: This is a workaround because the correct type of 
    !          R % Interactions is not being recognized in K % Initialize
    K % Interactions  =>  UNIVERSE % Interactions_NM_G

    end select !-- R
    end select !-- K
    end select !-- S

  end subroutine SetSlope_NM_G_I


  subroutine SetSlope_NM_G_DFV_I ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_NM_G_DFV_I_Form :: K )
    select type ( K )
      class is ( Slope_NM_G_DFV_I_Form )
    select type ( R  =>  S % CurrentSet )
      class is ( NeutrinoMoments_G_Form )

    call K % Initialize &
           ( S % RiemannSolver, S % DiffusionFactor, S % DivergenceTotal, R )
            !, IgnorabilityOption = S % IGNORABILITY )

    !-- FIXME: This is a workaround because the correct type of 
    !          R % Interactions is not being recognized in K % Initialize
    select type ( K2  =>  K % Component ( 2 ) % Element )
    class is ( Slope_NM_G_I_Form )
      K2 % Interactions  =>  UNIVERSE % Interactions_NM_G
      select type ( I  =>  UNIVERSE % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
        K2 % Communicator_X_1D  =>  I % Communicator_X_1D
      end select !-- I
    end select !-- K2
    
    end select !-- R
    end select !-- K
    end select !-- S

  end subroutine SetSlope_NM_G_DFV_I


end module Universe_R_CC__Form
