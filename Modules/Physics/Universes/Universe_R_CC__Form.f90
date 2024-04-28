module Universe_R_CC__Form

  !-- Universe_Radiation_CentralCore__Form

  use Basics
  use Mathematics
  use Gravitations
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
      Compute_dT_RI_CGS
    procedure, public, pass ( U ) :: &
      Compute_dT_RT_CGS
    procedure, public, pass ( U ) :: &
      Compute_dT_RK_F_CGS
    procedure, public, pass ( U ) :: &
      Compute_dT_RK_R_CGS
    procedure, public, pass ( U ) :: &
      Compute_dT_IS_F_CGS
    procedure, public, pass ( U ) :: &
      Compute_dT_IS_R_CGS
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
      SetSlope_F_P_SS, &      
      SetSlope_F_P_DFV_N, &
      SetSlope_F_P_DFV_N_SS, &
      SetSlope_NM_G_I, &
      SetSlope_NM_G_I_I, &
      SetSlope_NM_G_DFV_I

      private :: &
        Compute_dT_RI_CGS_Kernel, &
        Compute_dT_RT_CGS_Kernel, &
        Compute_dT_RK_F_CGS_Kernel, &
        Compute_dT_RK_R_CGS_Kernel

    interface
    
      module subroutine Compute_dT_RI_CGS_Kernel &
               ( dT_E, dT_N, ProperCell, Q, R, E, N, E_Eq, N_Eq, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT_E, dT_N
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Q, R, &
          E, N, &
          E_Eq, N_Eq
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_RI_CGS_Kernel

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

      module subroutine Compute_dT_RK_F_CGS_Kernel &
               ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, ProperCell, &
                 E_E, E_S_1, E_S_2, E_S_3, E_N, E, S_1, S_2, S_3, N, &
                 dT, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          E_E, E_S_1, E_S_2, E_S_3, E_N, &
            E,   S_1,   S_2,   S_3,   N
        real ( KDR ), intent ( in ) :: &
          dT
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_RK_F_CGS_Kernel

      module subroutine Compute_dT_RK_R_CGS_Kernel &
               ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, ProperCell, &
                 E_E, E_S_1, E_S_2, E_S_3, E_N, E, S_1, S_2, S_3, N, &
                 J_RD, SF_RD, N_RD, DI, dT, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          E_E, E_S_1, E_S_2, E_S_3, E_N, &
            E,   S_1,   S_2,   S_3,   N, &   
          J_RD, SF_RD, N_RD, &
          DI
        real ( KDR ), intent ( in ) :: &
          dT
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_RK_R_CGS_Kernel

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
      nStages
    logical ( KDL ) :: &
      ImplicitExplicit
    character ( LDL ) :: &
      RiemannSolverType

    ImplicitExplicit  =  .true.
    call PROGRAM_HEADER % GetParameter ( ImplicitExplicit, 'ImplicitExplicit' )

    if ( ImplicitExplicit ) then
      nStages  =  3  !-- IMEX
    else
      nStages  =  2
    end if
    call PROGRAM_HEADER % GetParameter ( nStages, 'nStages' )

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

      if ( ImplicitExplicit ) then
        !-- SetSlopeExplicit set to DFV by default in Step_RK_CS__Form
        S_R % SetSlopeImplicit  =>  SetSlope_NM_G_I_I
      else
        S_R % SetSlopeExplicit  =>  SetSlope_NM_G_DFV_I
      end if

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

      if ( ImplicitExplicit ) then
        S_F % SetSlopeExplicit  =>  SetSlope_F_P_DFV_N
        S_F % SetSlopeImplicit  =>  SetSlope_F_P_SS
      else
        S_F % SetSlopeExplicit  =>  SetSlope_F_P_DFV_N_SS
      end if

      !-- Combined step
      call S % Initialize &
             ( R, F, &
               ImplicitExplicitOption = ImplicitExplicit, &
               nStagesOption = nStages )

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
    allocate ( I % dT_Label ( 13 ) )

    I % dT_Label (  1 )  =  'GravitationAcceleration'
    I % dT_Label (  2 )  =  'FluidAdvection'
    I % dT_Label (  3 )  =  'RadiationStreaming'
    I % dT_Label (  4 )  =  'FluidEnergyError'
    I % dT_Label (  5 )  =  'FluidMomentum_1_Error'
    I % dT_Label (  6 )  =  'FluidMomentum_2_Error'
    I % dT_Label (  7 )  =  'FluidMomentum_3_Error'
    I % dT_Label (  8 )  =  'FluidNumberError'
    I % dT_Label (  9 )  =  'RadiationEnergyError'
    I % dT_Label ( 10 )  =  'RadiationMomentum_1_Error'
    I % dT_Label ( 11 )  =  'RadiationMomentum_2_Error'
    I % dT_Label ( 12 )  =  'RadiationMomentum_3_Error'
    I % dT_Label ( 13 )  =  'RadiationNumberError'

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


  subroutine Compute_dT_RI_CGS ( dT_E, dT_N, U, iC, T_Option )

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
    select type ( S_I  =>  S_R % SlopeSumExplicit % Component ( 2 ) % Element )
      class is ( Slope_NM_G_I_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( NeutrinoMoments_G_Form )
    select type ( A  =>  R % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (   C  =>  A   % Chart_GS, &
        SIV  =>  S_I % Storage_GS % Value, &
         RV  =>  R   % Storage_GS % Value )

    call Search ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_B )
    call Search ( R % iaBalanced, R % NUMBER_DENSITY_B, iNumber_B )

    call Compute_dT_RI_CGS_Kernel &
           ( dT_E, dT_N, C % ProperCell, &
             Q     =  SIV ( :, iEnergy_B ), &
             R     =  SIV ( :, iNumber_B ), &
             E     =   RV ( :, R % ENERGY_DENSITY_C ), &
             N     =   RV ( :, R % NUMBER_DENSITY_C ), &
             E_Eq  =   RV ( :, R % ENERGY_DENSITY_C_EQ ), &
             N_Eq  =   RV ( :, R % NUMBER_DENSITY_C_EQ ), &
             UseDeviceOption = R % DeviceMemory )

    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Universe_R_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_dT_RI_CGS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end select !-- R
    end select !-- S_R_I
    end select !-- S_R
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_RI_CGS


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
    select type ( S_F  =>  S % Step_CS_2 )
      class is ( Step_RK_CS_Form )
    select type ( S_F_SS  =>  S_F % SlopeSumExplicit % Component ( 2 ) % Element )
      class is ( Slope_F_P_SS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_HN_Form )
!    associate &
!      ( FS =>  F % SplitSource )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (   C  =>  A      % Chart_GS, &
        FSV  =>  S_F_SS % Storage_GS % Value, &
         FV  =>  F      % Storage_GS % Value )

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

 !   end associate !-- FS
    end select !-- F
    end select !-- S_F_SS
    end select !-- S_F
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_RT_CGS


  subroutine Compute_dT_RK_F_CGS &
               ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, U, iC, T_Option )

    real ( KDR ), intent ( inout ) :: &
      dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N
    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iEnergy_B, iNumber_B, &
      iMomentum_B_1, iMomentum_B_2, iMomentum_B_3

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_CS_Form )
    select type ( S_F  =>  S % Step_CS_2 )
      class is ( Step_RK_CS_Form )
    associate &
      ( E  =>  S_F % Error )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_HN_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (   C  =>  A % Chart_GS, &
         EV  =>  E % Storage_GS % Value, &
         FV  =>  F % Storage_GS % Value )

    call Search ( F % iaBalanced, F % ENERGY_DENSITY_B,     iEnergy_B )
    call Search ( F % iaBalanced, F % ELECTRON_DENSITY_B,   iNumber_B )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_B_1 )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_B_2 )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_B_3 )

    call Compute_dT_RK_F_CGS_Kernel &
           ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, C % ProperCell, &
             E_E    =  EV ( :, iEnergy_B ), &
             E_S_1  =  EV ( :, iMomentum_B_1 ), &
             E_S_2  =  EV ( :, iMomentum_B_2 ), &
             E_S_3  =  EV ( :, iMomentum_B_3 ), &
             E_N    =  EV ( :, iNumber_B ), &
             E      =  FV ( :, F % ENERGY_DENSITY_B ), &
             S_1    =  FV ( :, F % MOMENTUM_DENSITY_D_1 ), &
             S_2    =  FV ( :, F % MOMENTUM_DENSITY_D_2 ), &
             S_3    =  FV ( :, F % MOMENTUM_DENSITY_D_3 ), &
             N      =  FV ( :, F % ELECTRON_DENSITY_B ), &
             dT     =  I % dT  /  I % RampFactor, &
             UseDeviceOption = F % DeviceMemory )

    end associate !-- C, etc.
    end select !-- A
    end select !-- F
    end associate !-- E
    end select !-- S_F
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_RK_F_CGS


  subroutine Compute_dT_RK_R_CGS &
               ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, U, iC, T_Option )

    real ( KDR ), intent ( inout ) :: &
      dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N
    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iEnergy_B, iNumber_B, &
      iMomentum_B_1, iMomentum_B_2, iMomentum_B_3

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_CS_Form )
    select type ( S_R  =>  S % Step_CS_1 )
      class is ( Step_RK_CS_Form )
    associate &
      ( E  =>  S_R % Error )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( NeutrinoMoments_G_Form )
    select type ( A  =>  R % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      (   C  =>  A % Chart_GS, &
         EV  =>  E % Storage_GS % Value, &
         RV  =>  R % Storage_GS % Value )

    call Search ( R % iaBalanced, R % ENERGY_DENSITY_B,       iEnergy_B )
    call Search ( R % iaBalanced, R % NUMBER_DENSITY_B,       iNumber_B )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_B_1 )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_B_2 )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_B_3 )

    call Compute_dT_RK_R_CGS_Kernel &
           ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, C % ProperCell, &
             E_E    =  EV ( :, iEnergy_B ), &
             E_S_1  =  EV ( :, iMomentum_B_1 ), &
             E_S_2  =  EV ( :, iMomentum_B_2 ), &
             E_S_3  =  EV ( :, iMomentum_B_3 ), &
             E_N    =  EV ( :, iNumber_B ), &
             E      =  RV ( :, R % ENERGY_DENSITY_B ), &
             S_1    =  RV ( :, R % MOMENTUM_DENSITY_B_D_1 ), &
             S_2    =  RV ( :, R % MOMENTUM_DENSITY_B_D_2 ), &
             S_3    =  RV ( :, R % MOMENTUM_DENSITY_B_D_3 ), &
             N      =  RV ( :, R % NUMBER_DENSITY_B ), &
             J_RD   =  RV ( :, R % ENERGY_DENSITY_C_RD ), &
             SF_RD  =  RV ( :, R % STRESS_FACTOR_RD ), &
             N_RD   =  RV ( :, R % NUMBER_DENSITY_C_RD ), &
             DI     =  RV ( :, R % DIFFUSION_INDICATOR ), &
             dT     =  I % dT  /  I % RampFactor, &
             UseDeviceOption = R % DeviceMemory )

    end associate !-- C, etc.
    end select !-- A
    end select !-- R
    end associate !-- E
    end select !-- S_R
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_RK_R_CGS


  subroutine Compute_dT_IS_F_CGS &
               ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, U, dT_RS, iC, T_Option )

    real ( KDR ), intent ( inout ) :: &
      dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N
    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ) :: &
      iC
    real ( KDR ), intent ( in ) :: &
      dT_RS  !-- RadiationStreaming
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iEnergy_B, &
      iMomentum_B_1, iMomentum_B_2, iMomentum_B_3, &
      iNumber_B
    real ( KDR ) :: &
      FactorPoor, &
      FactorFair, &
      FactorGood, &
      FactorExcellent, &
      Factor

    FactorPoor       =  0.5_KDR
    FactorFair       =  0.8_KDR
    FactorGood       =  1.0_KDR
    FactorExcellent  =  1.1_KDR
    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_HN_Form )

    call Search ( F % iaBalanced, F % ENERGY_DENSITY_B,       iEnergy_B )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_B_1 )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_B_2 )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_B_3 )
    call Search ( F % iaBalanced, F % ELECTRON_DENSITY_B,   iNumber_B )

    associate ( IQ_E  =>  S % ImplicitQuality_2 ( iEnergy_B, : ) )
    Factor  =  1.0_KDR
    if ( any ( IQ_E  ==  S % IMPLICIT_POOR ) ) then
      Factor  =  FactorPoor
    else if ( any ( IQ_E  ==  S % IMPLICIT_FAIR ) ) then
      Factor  =  FactorFair
    else if ( any ( IQ_E  ==  S % IMPLICIT_GOOD ) ) then
      Factor  =  FactorGood
    else if ( all ( IQ_E  ==  S % IMPLICIT_EXCELLENT ) ) then
      Factor  =  FactorExcellent
    end if
    end associate !-- IQ_E
    dT_E  =  Factor  *  ( I % dT  /  I % RampFactor )
!    dT_E  =  max ( Factor  *  ( I % dT  /  I % RampFactor ), &
!                   0.1_KDR  *  dT_RS )

    ! associate ( IQ_S_1  =>  S % ImplicitQuality_2 ( iMomentum_B_1, : ) )
    ! Factor  =  1.0_KDR
    ! if ( any ( IQ_S_1  ==  S % IMPLICIT_POOR ) ) then
    !   Factor  =  FactorPoor
    ! ! else if ( any ( IQ_S_1  ==  S % IMPLICIT_FAIR ) ) then
    ! !   Factor  =  FactorFair
    ! ! else if ( any ( IQ_S_1  ==  S % IMPLICIT_GOOD ) ) then
    ! !   Factor  =  FactorGood
    ! else if ( all ( IQ_S_1  ==  S % IMPLICIT_EXCELLENT ) ) then
    !   Factor  =  FactorExcellent
    ! end if
    ! end associate !-- IQ_S_1
    ! dT_S_1  =  Factor  *  ( I % dT  /  I % RampFactor )

    ! associate ( IQ_S_2  =>  S % ImplicitQuality_2 ( iMomentum_B_2, : ) )
    ! Factor  =  1.0_KDR
    ! if ( any ( IQ_S_2  ==  S % IMPLICIT_POOR ) ) then
    !   Factor  =  FactorPoor
    ! ! else if ( any ( IQ_S_2  ==  S % IMPLICIT_FAIR ) ) then
    ! !   Factor  =  FactorFair
    ! ! else if ( any ( IQ_S_2  ==  S % IMPLICIT_GOOD ) ) then
    ! !   Factor  =  FactorGood
    ! else if ( all ( IQ_S_2  ==  S % IMPLICIT_EXCELLENT ) ) then
    !   Factor  =  FactorExcellent
    ! end if
    ! end associate !-- IQ_S_1
    ! dT_S_2  =  Factor  *  ( I % dT  /  I % RampFactor )

    ! associate ( IQ_S_3  =>  S % ImplicitQuality_2 ( iMomentum_B_3, : ) )
    ! Factor  =  1.0_KDR
    ! if ( any ( IQ_S_3  ==  S % IMPLICIT_POOR ) ) then
    !   Factor  =  FactorPoor
    ! ! else if ( any ( IQ_S_3  ==  S % IMPLICIT_FAIR ) ) then
    ! !   Factor  =  FactorFair
    ! ! else if ( any ( IQ_S_3  ==  S % IMPLICIT_GOOD ) ) then
    ! !   Factor  =  FactorGood
    ! else if ( all ( IQ_S_3  ==  S % IMPLICIT_EXCELLENT ) ) then
    !   Factor  =  FactorExcellent
    ! end if
    ! end associate !-- IQ_S_3
    ! dT_S_3  =  Factor  *  ( I % dT  /  I % RampFactor )

    associate ( IQ_N  =>  S % ImplicitQuality_2 ( iNumber_B, : ) )
    Factor  =  1.0_KDR
    if ( any ( IQ_N  ==  S % IMPLICIT_POOR ) ) then
      Factor  =  FactorPoor
    ! else if ( any ( IQ_N  ==  S % IMPLICIT_FAIR ) ) then
    !   Factor  =  FactorFair
    ! else if ( any ( IQ_N  ==  S % IMPLICIT_GOOD ) ) then
    !   Factor  =  FactorGood
    else if ( all ( IQ_N  ==  S % IMPLICIT_EXCELLENT ) ) then
      Factor  =  FactorExcellent
    end if
    end associate !-- IQ_N
    dT_N  =  Factor  *  ( I % dT  /  I % RampFactor )
!    dT_N  =  max ( Factor  *  ( I % dT  /  I % RampFactor ), &
!                   0.1_KDR  *  dT_RS )

    end select !-- F
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_IS_F_CGS


  subroutine Compute_dT_IS_R_CGS &
               ( dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N, U, iC, T_Option )

    real ( KDR ), intent ( inout ) :: &
      dT_E, dT_S_1, dT_S_2, dT_S_3, dT_N
    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iEnergy_B, &
      iMomentum_B_1, iMomentum_B_2, iMomentum_B_3, &
      iNumber_B
    real ( KDR ) :: &
      FactorPoor, &
      FactorFair, &
      FactorGood, &
      FactorExcellent, &
      Factor

    FactorPoor       =  0.5_KDR
    FactorFair       =  0.8_KDR
    FactorGood       =  1.0_KDR
    FactorExcellent  =  1.1_KDR
    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( NeutrinoMoments_G_Form )

    call Search ( R % iaBalanced, R % ENERGY_DENSITY_B,       iEnergy_B )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_B_1 )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_B_2 )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_B_3 )
    call Search ( R % iaBalanced, R % NUMBER_DENSITY_B,       iNumber_B )

    associate ( IQ_E  =>  S % ImplicitQuality_1 ( iEnergy_B, : ) )
    Factor  =  1.0_KDR
    if ( any ( IQ_E  ==  S % IMPLICIT_POOR ) ) then
      Factor  =  FactorPoor
    ! else if ( any ( IQ_E  ==  S % IMPLICIT_FAIR ) ) then
    !   Factor  =  FactorFair
    ! else if ( any ( IQ_E  ==  S % IMPLICIT_GOOD ) ) then
    !   Factor  =  FactorGood
    else if ( all ( IQ_E  ==  S % IMPLICIT_EXCELLENT ) ) then
      Factor  =  FactorExcellent
    end if
    end associate !-- IQ_E
    dT_E  =  Factor  *  ( I % dT  /  I % RampFactor )

    ! associate ( IQ_S_1  =>  S % ImplicitQuality_1 ( iMomentum_B_1, : ) )
    ! Factor  =  1.0_KDR
    ! if ( any ( IQ_S_1  ==  S % IMPLICIT_POOR ) ) then
    !   Factor  =  FactorPoor
    ! ! else if ( any ( IQ_S_1  ==  S % IMPLICIT_FAIR ) ) then
    ! !   Factor  =  FactorFair
    ! ! else if ( any ( IQ_S_1  ==  S % IMPLICIT_GOOD ) ) then
    ! !   Factor  =  FactorGood
    ! else if ( all ( IQ_S_1  ==  S % IMPLICIT_EXCELLENT ) ) then
    !   Factor  =  FactorExcellent
    ! end if
    ! end associate !-- IQ_S_1
    ! dT_S_1  =  Factor  *  ( I % dT  /  I % RampFactor )

    ! associate ( IQ_S_2  =>  S % ImplicitQuality_1 ( iMomentum_B_2, : ) )
    ! Factor  =  1.0_KDR
    ! if ( any ( IQ_S_2  ==  S % IMPLICIT_POOR ) ) then
    !   Factor  =  FactorPoor
    ! ! else if ( any ( IQ_S_2  ==  S % IMPLICIT_FAIR ) ) then
    ! !   Factor  =  FactorFair
    ! ! else if ( any ( IQ_S_2  ==  S % IMPLICIT_GOOD ) ) then
    ! !   Factor  =  FactorGood
    ! else if ( all ( IQ_S_2  ==  S % IMPLICIT_EXCELLENT ) ) then
    !   Factor  =  FactorExcellent
    ! end if
    ! end associate !-- IQ_S_1
    ! dT_S_2  =  Factor  *  ( I % dT  /  I % RampFactor )

    ! associate ( IQ_S_3  =>  S % ImplicitQuality_1 ( iMomentum_B_3, : ) )
    ! Factor  =  1.0_KDR
    ! if ( any ( IQ_S_3  ==  S % IMPLICIT_POOR ) ) then
    !   Factor  =  FactorPoor
    ! ! else if ( any ( IQ_S_3  ==  S % IMPLICIT_FAIR ) ) then
    ! !   Factor  =  FactorFair
    ! ! else if ( any ( IQ_S_3  ==  S % IMPLICIT_GOOD ) ) then
    ! !   Factor  =  FactorGood
    ! else if ( all ( IQ_S_3  ==  S % IMPLICIT_EXCELLENT ) ) then
    !   Factor  =  FactorExcellent
    ! end if
    ! end associate !-- IQ_S_3
    ! dT_S_3  =  Factor  *  ( I % dT  /  I % RampFactor )

    associate ( IQ_N  =>  S % ImplicitQuality_1 ( iNumber_B, : ) )
    Factor  =  1.0_KDR
    if ( any ( IQ_N  ==  S % IMPLICIT_POOR ) ) then
      Factor  =  FactorPoor
    ! else if ( any ( IQ_N  ==  S % IMPLICIT_FAIR ) ) then
    !   Factor  =  FactorFair
    ! else if ( any ( IQ_N  ==  S % IMPLICIT_GOOD ) ) then
    !   Factor  =  FactorGood
    else if ( all ( IQ_N  ==  S % IMPLICIT_EXCELLENT ) ) then
      Factor  =  FactorExcellent
    end if
    end associate !-- IQ_N
    dT_N  =  Factor  *  ( I % dT  /  I % RampFactor )

    end select !-- R
    end select !-- S
    end select !-- I

  end subroutine Compute_dT_IS_R_CGS


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
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_H_Form )
    associate &  !-- See InitializeIntegrator subroutine herein
      ( dT_1   =>  dT_Candidate (  1 ), &
        dT_2   =>  dT_Candidate (  2 ), &
        dT_3   =>  dT_Candidate (  3 ), &
        dT_4   =>  dT_Candidate (  4 ), &
        dT_5   =>  dT_Candidate (  5 ), &
        dT_6   =>  dT_Candidate (  6 ), &
        dT_7   =>  dT_Candidate (  7 ), &
        dT_8   =>  dT_Candidate (  8 ), &
        dT_9   =>  dT_Candidate (  9 ), &
        dT_10  =>  dT_Candidate ( 10 ), &
        dT_11  =>  dT_Candidate ( 11 ), &
        dT_12  =>  dT_Candidate ( 12 ), &
        dT_13  =>  dT_Candidate ( 13 ) )

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

    ! !-- Radiative transfer steps

    ! call U % Compute_dT_RT_CGS ( dT_4, dT_5, iC, T_Option )
    ! dT_4  =  U % InteractionFactor  *  dT_4
    ! dT_5  =  U % InteractionFactor  *  dT_5

!    !-- Radiation interaction steps

!    call U % Compute_dT_RI_CGS ( dT_6, dT_7, iC, T_Option )
!    dT_6  =  U % InteractionFactor  *  dT_6
!    dT_7  =  U % InteractionFactor  *  dT_7

    if ( S % EmbeddedMethod ) then

      !-- Fluid error steps

!      if ( I % iCheckpoint  >  1 ) &
        call U % Compute_dT_RK_F_CGS &
               ( dT_4, dT_5, dT_6, dT_7, dT_8, iC, T_Option )

      !-- Radiation error steps

!      if ( I % iCheckpoint  >  1 ) &
        call U % Compute_dT_RK_R_CGS &
               ( dT_9, dT_10, dT_11, dT_12, dT_13, iC, T_Option )

    else if ( S % ImplicitExplicit ) then

      !-- Fluid implicit solver steps

!      if ( I % iCheckpoint  >  1 ) &
        call U % Compute_dT_IS_F_CGS &
               ( dT_4, dT_5, dT_6, dT_7, dT_8, dT_3, iC, T_Option )

      !-- Radiation implicit solver steps

!      if ( I % iCheckpoint  >  1 ) &
        call U % Compute_dT_IS_R_CGS &
               ( dT_9, dT_10, dT_11, dT_12, dT_13, iC, T_Option )

    end if

   !-- Reduce across radiation types

   call CO % Initialize &
          ( I % Communicator_X_1D, nOutgoing = [ 6 ], &
            nIncoming = [ 6 ] )

   CO % Outgoing % Value ( 1 )      =  I % dT_Candidate ( 3 )
   CO % Outgoing % Value ( 2 : 6 )  =  I % dT_Candidate ( 9 : 13 )
   call CO % Reduce ( REDUCTION % MIN )
   I % dT_Candidate ( 3 )       =  CO % Incoming % Value ( 1 )
   I % dT_Candidate ( 9 : 13 )  =  CO % Incoming % Value ( 2 : 6 )

    end associate !-- dT_1, etc.
    end select !-- S
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


  subroutine SetSlope_F_P_SS ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_F_P_SS_Form :: K )
    select type ( K )
      class is ( Slope_F_P_SS_Form )
    select type ( F  =>  S % CurrentSet )
      class is ( Fluid_P_Form )

    call K % Initialize ( F )!, &
!             IgnorabilityOption = S % IGNORABILITY )

    end select !-- F
    end select !-- K
    end select !-- S

  end subroutine SetSlope_F_P_SS


  subroutine SetSlope_F_P_DFV_N ( S, K )

    !-- Compare SetSlope_N in Universe_F_C__Form

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

!    select type ( S )
!      class is ( Step_RK_CS_Form )
!    select type ( F  =>  S % CurrentSet )
!    class is ( Fluid_P_HN_Form )
!      allocate ( Slope_DFV_N_F_P_HN_Form :: K )
!    class default
      allocate ( Slope_DFV_N_Form :: K )
!    end select !-- F
!    end select !-- S

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
      call Search ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_B )
    end select !-- F

    if ( allocated ( S % DivergenceTotal ) ) then
      call K % Initialize &
             ( S % RiemannSolver, &
               S % DiffusionFactor, &
               S % DivergenceTotal, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B )
    else if ( allocated ( S % DivergencePart ) ) then
      call K % Initialize &
             ( S % RiemannSolver, &
               S % DiffusionFactor, &
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

  end subroutine SetSlope_F_P_DFV_N


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
    select type ( I  =>  UNIVERSE % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )
      K % Communicator_X_1D  =>  I % Communicator_X_1D
    end select !-- I

    end select !-- R
    end select !-- K
    end select !-- S

  end subroutine SetSlope_NM_G_I


  subroutine SetSlope_NM_G_I_I ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_NM_G_I_I_Form :: K )
    select type ( K )
      class is ( Slope_NM_G_I_I_Form )
    select type ( R  =>  S % CurrentSet )
      class is ( NeutrinoMoments_G_Form )

!     select type ( I  =>  R % Interactions )
!     class is ( Interactions_BM_Form )
! call Show ( '>>> Interactions_BM SetSlope_NM_G_I_I' )
!     end select

!     select type ( I  =>  R % Interactions )
!     type is ( Interactions_NM_G_Form )
! call Show ( '>>> Interactions_NM_G SetSlope_NM_G_I_I' )
!     end select

!     select type ( I  =>  UNIVERSE % Interactions_NM_G )
!     class is ( Interactions_NM_G_Form )
! call Show ( '>>> Interactions_NM_G SetSlope_NM_G_I_I UNIVERSE' )
!     end select

    call K % Initialize &
           ( R )!, &
!             IgnorabilityOption = S % IGNORABILITY )

    !-- FIXME: This is a workaround because the correct type of 
    !          R % Interactions is not being recognized in K % Initialize
    K % Interactions  =>  UNIVERSE % Interactions_NM_G
    select type ( I  =>  UNIVERSE % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )
      K % Communicator_X_1D  =>  I % Communicator_X_1D
    end select !-- I

    end select !-- R
    end select !-- K
    end select !-- S

  end subroutine SetSlope_NM_G_I_I


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
