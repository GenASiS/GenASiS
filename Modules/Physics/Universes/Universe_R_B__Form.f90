module Universe_R_B__Form

  !-- Universe_Radiation_Box__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Universe_F_B__Form

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: Universe_R_B_Form
    integer ( KDI ) :: &
      iRadiation  = 0, &
      nRadiations = 0
    real ( KDR ) :: &
      InteractionFactor
    logical ( KDL ) :: &
      ApplyStreaming, &
      ApplyInteractions, &
      EvolveFluid
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
    type ( CommunicatorForm ), allocatable :: &
      Communicator_PS  !-- PositionSpace
    type ( CollectiveOperation_R_Form ), dimension ( : ), allocatable :: &
      CO_SplitSource
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
    class ( Interactions_BM_Form ), allocatable :: &
      Interactions_BM
  contains
    procedure, private, pass :: &
      Initialize_R_B
    generic, public :: &
      Initialize => Initialize_R_B
    final :: &
      Finalize
    procedure, private, pass :: &
      SetCommunicator
    procedure, private, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeInteractions
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      InitializeSteps
    procedure, public, pass :: &
      InitializeIntegrator
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass ( U ) :: &
      Compute_dT_ET_CGS
  end type Universe_R_B_Form

    private :: &
      ResolveCycle_R, &
      PrepareStep_F, &
      Compute_dT_Local, &
      SetSlope_F_P_DFV_SS, &
      SetSlope_F_P_SS, &
      SetSlope_RM_I, &
      SetSlope_RM_DFV_I

      private :: &
        Compute_dT_ET_CGS_Kernel

    interface
    
      module subroutine Compute_dT_ET_CGS_Kernel &
               ( dT, ProperCell, Q, E, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Q, E
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_ET_CGS_Kernel

    end interface


contains


  subroutine Initialize_R_B &
               ( U, RadiationName, RadiationType, FormalismType, Name, &
                 ApplyStreamingOption, ApplyInteractionsOption, &
                 EvolveFluidOption, MinCoordinateOption, MaxCoordinateOption, &
                 FinishTimeOption, nCellsPositionOption, nWriteOption )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      ApplyStreamingOption, &
      ApplyInteractionsOption, &
      EvolveFluidOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsPositionOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_B'

    call U % Universe_H_Form % Initialize ( Name )

    !-- Radiations

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    !-- Operators

    U % ApplyStreaming    = .true.
    U % ApplyInteractions = .true.
    U % EvolveFluid       = .true.
    if ( present ( ApplyStreamingOption ) ) &
      U % ApplyStreaming = ApplyStreamingOption
    if ( present ( ApplyInteractionsOption ) ) &
      U % ApplyInteractions = ApplyInteractionsOption
    if ( present ( EvolveFluidOption ) ) &
      U % EvolveFluid = EvolveFluidOption

    !-- Units

    if ( .not. allocated ( U % Units_F ) ) then
      allocate ( U % Units_F ( 1 ) )
      call U % Units_F ( 1 ) % Initialize ( )
    end if

    if ( .not. allocated ( U % Units_R ) ) then
      allocate ( U % Units_R ( 1 ) )
      call U % Units_R ( 1 ) % Initialize ( )
    end if

    !-- Initializations

    call U % SetCommunicator &
           ( )
    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsPositionOption )
    ! call RB % InitializeMomentumSpace &
    !        ( EnergySpacingOption = EnergySpacingOption, &
    !          MinEnergyOption = MinEnergyOption, &
    !          MaxEnergyOption = MaxEnergyOption, &
    !          MinWidthEnergyOption = MinWidthEnergyOption, &
    !          EnergyScaleOption = EnergyScaleOption, &
    !          nCellsEnergyOption = nCellsEnergyOption )
    call U % InitializeGravitation &
           ( GravitationType = 'GALILEO' )
    call U % InitializeFluid &
           ( FluidType = 'IDEAL' )
    call U % InitializeInteractions &
           ( )
    call U % InitializeRadiation &
           ( )
    call U % InitializeSteps &
           ( )
    call U % InitializeIntegrator &
           ( FinishTimeOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % ResolveCycle      =>  ResolveCycle_R
    I % PrepareStep       =>  PrepareStep_F
    I % Compute_dT_Local  =>  Compute_dT_Local
    end associate !-- I

  end subroutine Initialize_R_B


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Interactions_BM ) ) &
      deallocate ( U % Interactions_BM )
    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
    if ( allocated ( U % CO_SplitSource ) ) &
      deallocate ( U % CO_SplitSource )
    if ( allocated ( U % Communicator_PS ) ) &
      deallocate ( U % Communicator_PS )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine SetCommunicator ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
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
      call Show ( 'Universe_R_B_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCommunicator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine SetCommunicator


  subroutine AllocateIntegrator ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )
      allocate ( Integrator_CS_1D_BM_CS_Form :: U % Integrator )
    case ( 'SPECTRAL' )
      allocate ( Integrator_CS_1D_CB_CS_Form :: U % Integrator )
    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_B_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine AllocateIntegrator


  subroutine InitializePositionSpace &
               ( U, CommunicatorOption, MinCoordinateOption, &
                 MaxCoordinateOption, nCellsOption )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption

    integer ( KDI ), dimension ( 3 ) :: &
      nCellsPosition

    nCellsPosition = [ 128, 128, 128 ]
    if ( present ( nCellsOption ) ) &
      nCellsPosition = nCellsOption
    call PROGRAM_HEADER % GetParameter ( nCellsPosition, 'nCellsPosition' )

    call U % Universe_F_B_Form % InitializePositionSpace &
           ( CommunicatorOption = U % Communicator_PS, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsPosition )

  end subroutine InitializePositionSpace


  subroutine InitializeInteractions ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )

    if ( allocated ( U % Interactions_BM ) ) &
      call U % Interactions_BM % Initialize ( F, U % Units_R )

    end select !-- F
    end select !-- I

  end subroutine InitializeInteractions


  subroutine InitializeRadiation ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )

      associate &
        ( G   =>  I % Geometry_X, &
          iR  =>  U % iRadiation )

      select case ( trim ( U % RadiationType ( iR ) ) )
      case ( 'GENERIC' )

        allocate ( RadiationMoments_BM_Form :: I % CurrentSet_X_1D )
        select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )

        call R % Initialize &
               ( G, U % Units_R, NameOption = U % RadiationName ( iR ) )
        if ( allocated ( U % Interactions_BM ) ) &
          call R % SetInteractions ( U % Interactions_BM )

        end select !-- R

      case ( 'PHOTONS' )

        allocate ( PhotonMoments_G_Form :: I % CurrentSet_X_1D )
        select type ( R  =>  I % CurrentSet_X_1D )
        class is ( PhotonMoments_G_Form )

        call R % Initialize &
               ( G, U % Units_R, NameOption = U % RadiationName ( iR ) )
        if ( allocated ( U % Interactions_BM ) ) &
          call R % SetInteractions ( U % Interactions_BM )

        end select !-- R

      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( U % RadiationType ( iR ), 'RadiationType', &
                    CONSOLE % ERROR )
        call Show ( 'Universe_R_B__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- RadiationType

      end associate !-- G, etc.

    end select !-- I

  end subroutine InitializeRadiation


  subroutine InitializeSteps ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    character ( LDL ) :: &
      RiemannSolverType

    !-- Fluid step

    if ( U % EvolveFluid ) then
 
      select type ( I  =>  U % Integrator )
        class is ( Integrator_CS_Form )
      associate &
        ( F  =>  I % CurrentSet_X )

      allocate ( Step_RK_CS_Form :: I % Step_X )
      select type ( S  =>  I % Step_X )
        class is ( Step_RK_CS_Form )

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

      S % SetSlope  =>  SetSlope_F_P_DFV_SS

      call S % Initialize ( F )

      end select !-- S
      end associate !-- F
      end select !-- I

    end if

    !-- Radiation step

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )

      associate &
        ( R  =>  I % CurrentSet_X_1D )

      allocate ( Step_RK_CS_Form :: I % Step_1D )
      select type ( S  =>  I % Step_1D )
        class is ( Step_RK_CS_Form )

      if ( U % ApplyStreaming .and. .not. U % ApplyInteractions ) then

        allocate ( DivergencePart_RM_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
        call DT % Initialize ( R )
        end associate !-- DT

      else if ( U % ApplyInteractions .and. .not. U % ApplyStreaming ) then

        S % SetSlope  =>  SetSlope_RM_I 

      else if ( U % ApplyStreaming .and. U % ApplyInteractions ) then

        allocate ( DivergencePart_RM_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
        call DT % Initialize ( R )
        end associate !-- DT

        allocate ( DiffusionFactor_RM_Form :: S % DiffusionFactor )
        associate ( DF  =>  S % DiffusionFactor )
        call DF % Initialize ( R )
        end associate !-- DF

        S % SetSlope  =>  SetSlope_RM_DFV_I

      end if !-- Radiation operators

      call S % Initialize ( R )

      end select !-- S
      end associate !-- R

    end select !-- I

  end subroutine InitializeSteps


  subroutine InitializeIntegrator ( U, FinishTimeOption, nWriteOption )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iCS  !-- iCurrentSet

    select type ( I => U % Integrator )
    class is ( Integrator_CS_1D_CS_Form )

    I % iCurrentSet  =  U % iRadiation

    I % nCurrentSets  =  size ( U % RadiationName )
    allocate ( I % dT_Label ( 3 ) )

    I % dT_Label ( 1 )  =  'FluidAdvection'
    I % dT_Label ( 2 )  =  'RadiationStreaming'
    I % dT_Label ( 3 )  =  'EnergyTransfer'

    U % InteractionFactor  =  1.0e-2_KDR  /  I % nCurrentSets
    call PROGRAM_HEADER % GetParameter &
           ( U % InteractionFactor, 'InteractionFactor' )

    I % StreamSuffix  =  '_' // trim ( U % RadiationName ( U % iRadiation ) )

    call I % Initialize &
           ( CommunicatorOption = U % Communicator_PS, &
             Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    end select !-- I

  end subroutine InitializeIntegrator


  subroutine ShowParameters ( U )

    class ( Universe_R_B_Form ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % RadiationName,     'RadiationName',     U % IGNORABILITY )
    call Show ( U % RadiationType,     'RadiationType',     U % IGNORABILITY )
    call Show ( U % iRadiation,        'iRadiation',        U % IGNORABILITY )
    call Show ( U % FormalismType,     'FormalismType',     U % IGNORABILITY )
    call Show ( U % ApplyStreaming,    'ApplyStreaming',    U % IGNORABILITY )
    call Show ( U % ApplyInteractions, 'ApplyInteractions', U % IGNORABILITY )
    call Show ( U % EvolveFluid,       'EvolveFluid',       U % IGNORABILITY )
    call Show ( U % InteractionFactor, 'InteractionFactor', U % IGNORABILITY )

  end subroutine ShowParameters


  subroutine Compute_dT_ET_CGS ( dT, U, iC, T_Option )

    real ( KDR ), intent ( inout ) :: &
      dT
    class ( Universe_R_B_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
  
    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( RadiationMoments_BM_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        FV  =>  F % Storage_GS % Value, &
        RV  =>  R % Storage_GS % Value )

    call Compute_dT_ET_CGS_Kernel &
           ( dT, C % ProperCell, &
             Q  =  RV ( :, R % HEATING_RATE ), &
             E  =  FV ( :, F % ENERGY_DENSITY_C ), &
             UseDeviceOption = F % DeviceMemory )

    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Universe_R_B_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_dT_ET_CGS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end select !-- F
    end select !-- R
    end select !-- I

  end subroutine Compute_dT_ET_CGS


  subroutine ResolveCycle_R ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    
    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( PhotonMoments_G_Form )

    call R % ComputeSpectralParameters ( )
    call R % ComputeHeatingRate ( )

    end select !-- R
    end select !-- I

  end subroutine ResolveCycle_R


  subroutine PrepareStep_F ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iEnergy_R, iEnergy_F, &
      nSources, &
      nValues
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, iMomentum_F
    real ( KDR ), dimension ( :, : ), pointer :: &
      RSB, &  !-- 2D alias for outgoing buffer
      FSB     !-- 2D alias for incoming buffer

    select type ( U  =>  I % System )
      class is ( Universe_R_B_Form )
    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( S_1D  =>  I % Step_1D )
      class is ( Step_RK_CS_Form )
    select type ( SS_R_I  =>  S_1D % SlopeSum % Component ( 2 ) % Element )
      class is ( Slope_RM_I_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( RadiationMoments_BM_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )

    call Search &
           ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )

    call Search &
           ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )

    nSources  =  4

    if ( .not. allocated ( U % CO_SplitSource ) ) &
      allocate ( U % CO_SplitSource ( F % Atlas % nCharts ) )

    do iC  =  1,  F % Atlas % nCharts
      associate &
        ( CO  =>  U % CO_SplitSource ( iC ) )
      associate &
        ( RSV  =>  SS_R_I % Storage ( iC ) % Value, &
          FSV  =>  F % SplitSource % Storage ( iC ) % Value )
      associate &
        ( RS_E    =>  RSV ( :, iEnergy_R ), &
          RS_S_1  =>  RSV ( :, iMomentum_R ( 1 ) ), &
          RS_S_2  =>  RSV ( :, iMomentum_R ( 2 ) ), &
          RS_S_3  =>  RSV ( :, iMomentum_R ( 3 ) ), &
          FS_G    =>  FSV ( :, iEnergy_F ), &
          FS_S_1  =>  FSV ( :, iMomentum_F ( 1 ) ), &
          FS_S_2  =>  FSV ( :, iMomentum_F ( 2 ) ), &
          FS_S_3  =>  FSV ( :, iMomentum_F ( 3 ) ) )

      nValues  =  size ( FSV, dim = 1 )

      if ( .not. allocated ( CO % Outgoing ) )  &
        call CO % Initialize &
               ( I % Communicator_X_1D, &
                 nOutgoing  =  [ nValues * nSources ], &
                 nIncoming  =  [ nValues * nSources ] )

      RSB ( 1 : nValues,  1 : nSources )  =>  CO % Outgoing % Value
      FSB ( 1 : nValues,  1 : nSources )  =>  CO % Incoming % Value

      call Copy ( RS_E,   RSB ( :, 1 ) )
      call Copy ( RS_S_1, RSB ( :, 2 ) )
      call Copy ( RS_S_2, RSB ( :, 3 ) )
      call Copy ( RS_S_3, RSB ( :, 4 ) )
      call Multiply ( RSB, -1.0_KDR )
      call CO % Reduce ( REDUCTION % SUM )
      call Copy ( FSB ( :, 1 ), FS_G )
      call Copy ( FSB ( :, 2 ), FS_S_1 )
      call Copy ( FSB ( :, 3 ), FS_S_2 )
      call Copy ( FSB ( :, 4 ), FS_S_3 )

      end associate !-- RS_E, etc.
      end associate !-- RSV, etc.
      end associate !-- CO
    end do !-- iC

    end select !-- F
    end select !-- R
    end select !-- SS_R_I
    end select !-- S_1D
    end select !-- I
    end select !-- U

  end subroutine PrepareStep_F


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

    associate &
      ( dT_1  =>  dT_Candidate ( 1 ), &
        dT_2  =>  dT_Candidate ( 2 ), &
        dT_3  =>  dT_Candidate ( 3 ) )

    select type ( U  =>  I % System )
      class is ( Universe_R_B_Form )
    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
     
      !-- CS ( Fluid )

      if ( U % EvolveFluid ) then
        call I % Compute_dT_CS_CGS &
               ( I % EigenspeedSet_X, dT_1, iC, T_Option )
        dT_1  =  I % CourantFactor  *  dT_1
      end if

      !-- CS_1D ( Radiation )

      if ( U % ApplyStreaming ) then
        call I % Compute_dT_CS_CGS &
               ( I % EigenspeedSet_X_1D, dT_2, iC, T_Option )
        dT_2  =  I % CourantFactor_1D  *  dT_2
      end if

      if ( U % ApplyInteractions ) then
        call U % Compute_dT_ET_CGS ( dT_3, iC, T_Option )
        dT_3  =  U % InteractionFactor  *  dT_3
      end if

      !-- Reduce across CS_1D

      call CO % Initialize &
             ( I % Communicator_X_1D, nOutgoing = [ 2 ], &
               nIncoming = [ 2 ] )

      CO % Outgoing % Value  =  I % dT_Candidate ( 2 : 3 )
      call CO % Reduce ( REDUCTION % MIN )
      I % dT_Candidate ( 2 : 3 )  =  CO % Incoming % Value

    end select !-- I
    end select !-- U
    end associate !-- dT_1, etc.

  end subroutine Compute_dT_Local


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


  subroutine SetSlope_F_P_DFV_SS ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_F_P_DFV_SS_Form :: K )
    select type ( K )
      class is ( Slope_F_P_DFV_SS_Form )
    select type ( F  =>  S % CurrentSet )
      class is ( Fluid_P_Form )

    call K % Initialize &
           ( S % RiemannSolver, S % DiffusionFactor, S % DivergenceTotal, F )
!             IgnorabilityOption = S % IGNORABILITY )

    end select !-- F
    end select !-- K
    end select !-- S

  end subroutine SetSlope_F_P_DFV_SS


  subroutine SetSlope_RM_I ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_RM_I_Form :: K )
    select type ( K )
      class is ( Slope_RM_I_Form )
    select type ( R  =>  S % CurrentSet )
      class is ( RadiationMoments_BM_Form )

    call K % Initialize &
           ( R )!, &
!             IgnorabilityOption = S % IGNORABILITY )

    end select !-- R
    end select !-- K
    end select !-- S

  end subroutine SetSlope_RM_I


  subroutine SetSlope_RM_DFV_I ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )

    allocate ( Slope_RM_DFV_I_Form :: K )
    select type ( K )
      class is ( Slope_RM_DFV_I_Form )
    select type ( R  =>  S % CurrentSet )
      class is ( RadiationMoments_BM_Form )

    call K % Initialize &
           ( S % RiemannSolver, S % DiffusionFactor, S % DivergenceTotal, R )
            !, IgnorabilityOption = S % IGNORABILITY )

    end select !-- R
    end select !-- K
    end select !-- S

  end subroutine SetSlope_RM_DFV_I


end module Universe_R_B__Form
