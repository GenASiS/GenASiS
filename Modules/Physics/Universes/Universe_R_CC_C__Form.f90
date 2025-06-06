module Universe_R_CC_C__Form

  !-- Universe_Radiation_CentralCore_Collected__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Radiations
  use Measures_R_CC_C__Form
  use Series_R_CC_C__Form
  use Universe_F_CC__Form

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: Universe_R_CC_C_Form
    integer ( KDI ) :: &
      nRadiations = 0
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
    type ( SphericalAverageForm ), dimension ( : ), allocatable :: &
      SA_Radiation, &
      SA_Interactions
    type ( AzimuthalAverageForm ), dimension ( : ), allocatable :: &
      AA_Radiation, &
      AA_Interactions
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
    type ( Coarsening_C_RM_Form ), dimension ( : ), allocatable :: &
      Coarsening_R
    class ( Interactions_NM_G_Form ), dimension ( : ), allocatable :: &
      Interactions_NM_G
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    final :: &
      Finalize
    procedure, public, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      InitializeInteractions
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, public, pass :: &
      InitializeStep
    procedure, public, pass :: &
      InitializeIntegrator
    procedure, public, pass :: &
      InitializeDiagnostics
    procedure, public, pass :: &
      SetMeasures
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass :: &
      ShowDiagnostics
  end type Universe_R_CC_C_Form

    private :: &
      Compute_dT_Local, &
      InitializeSeries, &
      Analyze, &
      Write, &
      Average, &
      Set_T_CheckpointInterval, &
      SetSlope_F_P_DFV_N


contains


  subroutine Initialize_R_CC_C &
               ( U, RadiationName, RadiationType, FormalismType, FluidType, &
                 GravitationType, Name, UnitsTypeOption, FinishTimeOption, &
                 RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_R_CC_C_Form ), intent ( inout ), target :: &
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
      U % Type  =  'a Universe_R_CC_C'

    call U % Universe_H_Form % Initialize &
           ( Name, UnitsTypeOption = UnitsTypeOption )

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
    
    !-- Radiations

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    !-- Initializations

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )
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
    call U % InitializeStep &
           ( )
    call U % InitializeIntegrator &
           ( GravitationType, &
             FinishTimeOption = FinishTimeOption, &
             nWriteOption = nWriteOption )
    call U % InitializeDiagnostics &
           ( )

    call U % SetMeasures ( )

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % Compute_dT_Local          =>  Compute_dT_Local
    I % InitializeSeries          =>  InitializeSeries
    I % Analyze                   =>  Analyze
    I % Write                     =>  Write
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    end associate !-- I

  end subroutine Initialize_R_CC_C


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Interactions_NM_G ) ) &
      deallocate ( U % Interactions_NM_G )
    if ( allocated ( U % Coarsening_R ) ) &
      deallocate ( U % Coarsening_R )
    if ( allocated ( U % AA_Interactions ) ) &
      deallocate ( U % AA_Interactions )
    if ( allocated ( U % AA_Radiation ) ) &
      deallocate ( U % AA_Radiation )
    if ( allocated ( U % SA_Interactions ) ) &
      deallocate ( U % SA_Interactions )
    if ( allocated ( U % SA_Radiation ) ) &
      deallocate ( U % SA_Radiation )
    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine AllocateIntegrator ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )
      allocate ( Integrator_CS_1D_C_CS_Form :: U % Integrator )
    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_CC_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine AllocateIntegrator


  subroutine InitializeRadiation ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

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

      associate ( nR  =>  U % nRadiations )

      allocate ( NeutrinoMoments_G_Form :: I % CurrentSet_X_1D ( nR ) )

      if ( allocated ( U % PositionSpace_AA ) ) &
        allocate ( U % AA_Radiation ( nR ) )
      if ( allocated ( U % PositionSpace_SA ) ) &
        allocate ( U % SA_Radiation ( nR ) )
    
      do iR  =  1, nR

        select case ( trim ( U % RadiationType ( iR ) ) )
        case ( 'NEUTRINOS_E', 'NEUTRINOS_EB', 'NEUTRINOS_HL' )

        select type ( R  =>  I % CurrentSet_X_1D ( iR ) )
        class is ( NeutrinoMoments_G_Form )

        call R % Initialize &
               ( F, U % Units_R, U % RadiationType ( iR ), &
                 NameOption = U % RadiationName ( iR ) )

        !-- Azimuthal average
        if ( allocated ( U % PositionSpace_AA ) ) then
          associate &
            ( AA     =>  U % AA_Radiation ( iR ), &
               A_AA  =>  U % PositionSpace_AA )
          allocate ( NeutrinoMoments_G_Form :: AA % FieldSet_AA )
          select type ( R_AA  =>  AA % FieldSet_AA )
            type is ( NeutrinoMoments_G_Form )
          select type ( F_AA  =>  U % AA_Fluid % FieldSet_AA )
            class is ( Fluid_P_HN_Form )
          call R_AA % Initialize &
                 ( F_AA, U % Units_R, R % RadiationType, &
                   NameOption = trim ( R % Name ) // '_AA' )
          call AA % Initialize &
                 ( G, R, A_AA, &
                   iaAverageOption = [ R % iaBalanced ] )
          end select !-- F_AA
          end select !-- R_AA
          end associate !-- AA, etc.
        end if !-- allocated PositionSpace_AA

        !-- Spherical average
        if ( allocated ( U % PositionSpace_SA ) ) then
          associate &
            ( SA     =>  U % SA_Radiation ( iR ), &
               A_SA  =>  U % PositionSpace_SA )
          allocate ( NeutrinoMoments_G_Form :: SA % FieldSet_SA )
          select type ( R_SA  =>  SA % FieldSet_SA )
            type is ( NeutrinoMoments_G_Form )
          select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
            class is ( Fluid_P_HN_Form )
          call R_SA % Initialize &
                 ( F_SA, U % Units_R, R % RadiationType, &
                   NameOption = trim ( R % Name ) // '_SA' )
          call SA % Initialize &
                 ( G, R, A_SA, &
                   iaAverageOption = [ R % iaBalanced ] )
          end select !-- F_SA
          end select !-- R_SA
          end associate !-- SA, etc.
        end if !-- allocated PositionSpace_SA

        end select !-- R

        case default
          call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
          call Show ( U % RadiationType ( iR ), 'RadiationType', &
                      CONSOLE % ERROR )
          call Show ( 'Universe_R_CC_C__Form', 'module', CONSOLE % ERROR )
          call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- RadiationType

      end do !-- iR

      end associate !-- nR
      end select !-- F

    end associate !-- G
    end select !-- I

  end subroutine InitializeRadiation


  subroutine InitializeInteractions ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR, &
      iRB
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( NeutrinoMoments_G_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )
    associate &
      (  G  =>  I % Geometry_X, &
        nR  =>  U % nRadiations )

    allocate ( iaAverage ( 0 ) )

    if ( allocated ( U % Interactions_NM_G ) ) then
      associate ( Int  =>  U % Interactions_NM_G )
      if ( allocated ( U % PositionSpace_AA ) ) &
        allocate ( U % AA_Interactions ( nR ) )
      if ( allocated ( U % PositionSpace_SA ) ) &
        allocate ( U % SA_Interactions ( nR ) )
      do iR  =  1,  nR

        select case ( iR )
        case ( 1 )  !-- NEUTRINOS_E
          iRB  =  2
        case ( 2 )  !-- NEUTRINOS_EB
          iRB  =  1
        case ( 3 )  !-- NEUTRINOS_HL
          iRB  =  3
        end select !-- iR

        call Int ( iR ) % Initialize &
               ( R ( iR ), R ( iRB ), U % Units_R, F )
        call R ( iR ) % SetInteractions ( Int ( iR ) )

        !-- Azimuthal average
        if ( allocated ( U % PositionSpace_AA ) ) then
          associate &
            ( AA     =>  U % AA_Interactions ( iR ), &
               A_AA  =>  U % PositionSpace_AA )
          allocate ( Interactions_NM_G_Form :: AA % FieldSet_AA )
          select type ( I_AA  =>  AA % FieldSet_AA )
            type is ( Interactions_NM_G_Form )
          select type ( R_AA  =>  U % AA_Radiation ( iR ) % FieldSet_AA )
            type is ( NeutrinoMoments_G_Form )
          select type ( RB_AA  =>  U % AA_Radiation ( iRB ) % FieldSet_AA )
            type is ( NeutrinoMoments_G_Form )
          select type ( F_AA  =>  U % AA_Fluid % FieldSet_AA )
            type is ( Fluid_P_HN_Form )
          call I_AA % Initialize &
                 ( R_AA, RB_AA, U % Units_R, F_AA, &
                   NameOption = trim ( Int ( iR ) % Name ) // '_AA' )
          call AA % Initialize &
                 ( G, Int ( iR ), A_AA, iaAverageOption = iaAverage )
          call R_AA % SetInteractions ( I_AA )
          end select !-- F_AA
          end select !-- RB_AA
          end select !-- R_AA
          end select !-- I_AA
          end associate !-- AA, etc.
        end if !-- allocated PositionSpace_AA

        !-- Spherical average
        if ( allocated ( U % PositionSpace_SA ) ) then
          associate &
            ( SA     =>  U % SA_Interactions ( iR ), &
               A_SA  =>  U % PositionSpace_SA )
          allocate ( Interactions_NM_G_Form :: SA % FieldSet_SA )
          select type ( I_SA  =>  SA % FieldSet_SA )
            type is ( Interactions_NM_G_Form )
          select type ( R_SA  =>  U % SA_Radiation ( iR ) % FieldSet_SA )
            type is ( NeutrinoMoments_G_Form )
          select type ( RB_SA  =>  U % SA_Radiation ( iRB ) % FieldSet_SA )
            type is ( NeutrinoMoments_G_Form )
          select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
            type is ( Fluid_P_HN_Form )
          call I_SA % Initialize &
                 ( R_SA, RB_SA, U % Units_R, F_SA, &
                   NameOption = trim ( Int ( iR ) % Name ) // '_SA' )
          call SA % Initialize &
                 ( G, Int ( iR ), A_SA, iaAverageOption = iaAverage )
          call R_SA % SetInteractions ( I_SA )
          end select !-- F_SA
          end select !-- RB_SA
          end select !-- R_SA
          end select !-- I_SA
          end associate !-- SA, etc.
        end if !-- allocated PositionSpace_SA

      end do !-- iR
      end associate !-- Int
    end if

    end associate !-- G, etc.
    end select !-- F
    end select !-- R
    end select !-- I

  end subroutine InitializeInteractions


  subroutine SetBoundaryConditions ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U
    
    integer ( KDI ) :: &
      iR

    call U % Universe_F_CC_Form % SetBoundaryConditions ( )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )

    associate &
      ( R  =>  I % CurrentSet_X_1D )

    do iR  =  1, U % nRadiations
      call R ( iR ) % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
      call R ( iR ) % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
      call R ( iR ) % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end do !-- iR

    end associate !-- R

    end select !-- I

  end subroutine SetBoundaryConditions


  subroutine InitializeStep ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR
    character ( LDL ) :: &
      RiemannSolverType

    U % Coarsen  =  .true.
    call PROGRAM_HEADER % GetParameter ( U % Coarsen, 'Coarsen' )

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_C_CS_Form )

      allocate ( Step_RK_NM_G_1D_C_Form :: I % Step_X )
      select type ( S  =>  I % Step_X )
        class is ( Step_RK_NM_G_1D_C_Form )

      associate ( nR  =>  U % nRadiations )

      allocate ( S % Step_CS_1D ( nR ) )
      allocate ( S % Step_CS )
      associate &
        ( S_R  =>  S % Step_CS_1D ( : ), &
          S_F  =>  S % Step_CS )

      !-- Radiation
      select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )

      do iR  =  1, nR
 
        allocate ( DivergencePart_NM_G_Form :: S_R ( iR ) % DivergenceTotal )
        associate ( DT  =>  S_R ( iR ) % DivergenceTotal )
        call DT % Initialize ( R ( iR ) )
        end associate !-- DT

        allocate ( DiffusionFactor_RM_Form :: S_R ( iR ) % DiffusionFactor )
        select type ( DF  =>  S_R ( iR ) % DiffusionFactor )
        class is ( DiffusionFactor_RM_Form )
          call DF % Initialize ( U % Interactions_NM_G ( iR ) )
        end select !-- DF

      end do !-- iR

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

      S_F % SetSlopeExplicit  =>  SetSlope_F_P_DFV_N

      !-- Combined step
      call S % Initialize &
             ( R, F, ImplicitExplicitOption = .true., nStagesOption = 3 )

      !-- Coarsening
      if ( U % Coarsen ) then

        allocate ( U % Coarsening_F )
        associate &
          ( C_F  =>  U % Coarsening_F, &
            G    =>  I % Geometry_X )
          call C_F % Initialize ( F, G )
          call S_F % SetCoarsening ( C_F )
        end associate !-- C_F, etc.

        allocate ( U % Coarsening_R ( nR ) )
        associate &
          ( C_R  =>  U % Coarsening_R, &
            G    =>  I % Geometry_X )
        do iR  =  1, nR
          call C_R ( iR ) % Initialize ( R ( iR ), G )
          call S_R ( iR ) % SetCoarsening ( C_R ( iR ) )
        end do !-- iR
        end associate !-- C_R, etc.
      end if

      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_R_CC_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select    !-- F
      end select    !-- R
      end associate !-- S_R, S_F
      end associate !-- nR
      end select    !-- S

    end select !-- I

  end subroutine InitializeStep


  subroutine InitializeIntegrator &
               ( U, GravitationType, FinishTimeOption, GravityFactorOption, &
                 nWriteOption )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      GravityFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iR

    select type ( I => U % Integrator )
    class is ( Integrator_CS_1D_C_CS_Form )

    allocate ( I % dT_Label ( 3 ) )

    I % dT_Label ( 1 )  =  'GravitationAcceleration'
    I % dT_Label ( 2 )  =  'FluidAdvection'
    I % dT_Label ( 3 )  =  'RadiationStreaming'

    U % GravityFactor  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter &
           ( U % GravityFactor, 'GravityFactor' )

    call I % Initialize &
           ( U % nRadiations, &
             Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    do iR  =  1, U % nRadiations
      call U % Interactions_NM_G ( iR ) % SetStream ( I % Checkpoint_X )
      call U % Interactions_NM_G ( iR ) % Show ( )
    end do

    end select !-- I

  end subroutine InitializeIntegrator


  subroutine InitializeDiagnostics ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    call U % Universe_F_CC_Form % InitializeDiagnostics ( )

    !-- AzimuthalAverage Stream

    if ( allocated ( U % PositionSpace_AA ) ) then
      associate ( S_AA  =>  U % Stream_AA )
      do iR  =  1,  U % nRadiations
        select type ( R_AA  =>  U % AA_Radiation ( iR ) % FieldSet_AA )
          class is ( NeutrinoMoments_G_Form )
        select type ( I_AA  =>  U % AA_Interactions ( iR ) % FieldSet_AA )
          class is ( Interactions_NM_G_Form )
        call R_AA % SetStream ( S_AA )
        call I_AA % SetStream ( S_AA )
        end select !-- I_AA
        end select !-- R_AA
      end do
      end associate !-- S_AA
    end if !-- allocated PositionSpace_AA

    !-- SphericalAverage Stream

    if ( allocated ( U % PositionSpace_SA ) ) then
      associate ( S_SA  =>  U % Stream_SA )
      do iR  =  1,  U % nRadiations
        select type ( R_SA  =>  U % SA_Radiation ( iR ) % FieldSet_SA )
          class is ( NeutrinoMoments_G_Form )
        select type ( I_SA  =>  U % SA_Interactions ( iR ) % FieldSet_SA )
          class is ( Interactions_NM_G_Form )
        call R_SA % SetStream ( S_SA )
        call I_SA % SetStream ( S_SA )
        end select !-- I_SA
        end select !-- R_SA
      end do
      end associate !-- S_SA
    end if !-- allocated PositionSpace_SA

  end subroutine InitializeDiagnostics


  subroutine SetMeasures ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ), target :: &
      U

    integer ( KDI ) :: &
      iR
    class ( Atlas_H_Form ), pointer :: &
      A_SA
    class ( FieldSet_BM_Form ), pointer :: &
      G_SA, &
      F_SA, &
      F
    type ( FieldSet_BM_Pointer ), dimension ( : ), allocatable :: &
      R_SA_1D

    associate ( nR  =>  U % nRadiations )
    allocate ( R_SA_1D ( nR ) )

    if ( allocated ( U % PositionSpace_SA ) ) then
      A_SA     =>  U % PositionSpace_SA
      G_SA     =>  U % SA_Gravitation % FieldSet_SA
      F_SA     =>  U % SA_Fluid % FieldSet_SA
      F        =>  U % SA_Fluid % FieldSet
      do iR  =  1, nR
        R_SA_1D ( iR ) % Pointer  =>  U % SA_Radiation ( iR ) % FieldSet_SA
      end do !-- iR
    else !-- 1D
      select type ( I  =>  U % Integrator )
        class is ( Integrator_CS_1D_C_CS_Form )
      A_SA     =>  I % X
      G_SA     =>  I % Geometry_X
      F_SA     =>  I % CurrentSet_X
      F        =>  I % CurrentSet_X
      do iR  =  1, nR
        R_SA_1D ( iR ) % Pointer  =>  I % CurrentSet_X_1D ( iR )
      end do !-- iR
      end select !-- I
    end if

    allocate ( Measures_R_CC_C_Form :: U % Measures )
    select type ( M  =>  U % Measures )
      class is ( Measures_R_CC_C_Form )
    call M % Initialize &
          ( R_SA_1D, F, F_SA, G_SA, A_SA, &
            Units_R = U % Units_R ( 1 ), Units_F = U % Units_F ( 1 ) )
    end select !-- M

    end associate !-- nR

  end subroutine SetMeasures


  subroutine ShowParameters ( U )

    class ( Universe_R_CC_C_Form ), intent ( in ) :: &
      U

    call U % Universe_F_CC_Form % ShowParameters ( )
    
    if ( allocated ( U % RadiationName ) ) &
      call Show ( U % RadiationName, 'RadiationName', U % IGNORABILITY )
    if ( allocated ( U % RadiationType ) ) &
      call Show ( U % RadiationType, 'RadiationType', U % IGNORABILITY )
    call Show ( U % FormalismType, 'FormalismType', U % IGNORABILITY )

  end subroutine ShowParameters


  subroutine ShowDiagnostics ( U )

      class ( Universe_R_CC_C_Form ), intent ( in ) :: &
        U

    integer ( KDI ) :: &
      iR

    if ( allocated ( U % Coarsening_F ) ) &
      call U % Coarsening_F % Show ( )
    do iR  =  1,  U % nRadiations
      if ( allocated ( U % Coarsening_R ) ) &
        call U % Coarsening_R ( iR ) % Show ( )
    end do !-- iR

    if ( allocated ( U % PositionSpace_AA ) ) then
      call U % PositionSpace_AA % Show ( )
      call U % AA_Gravitation % FieldSet_AA % Show ( )
      call U % AA_Fluid % FieldSet_AA % Show ( )
      do iR  =  1,  U % nRadiations
        call U % AA_Radiation ( iR ) % FieldSet_AA % Show ( )
        call U % AA_Interactions ( iR ) % FieldSet_AA % Show ( )
      end do !-- iR
      call U % Stream_AA % Show ( )
    end if !-- allocated PositionSpace_AA

    if ( allocated ( U % PositionSpace_SA ) ) then
      call U % PositionSpace_SA % Show ( )
      call U % SA_Gravitation % FieldSet_SA % Show ( )
      call U % SA_Fluid % FieldSet_SA % Show ( )
      do iR  =  1,  U % nRadiations
        call U % SA_Radiation ( iR ) % FieldSet_SA % Show ( )
        call U % SA_Interactions ( iR ) % FieldSet_SA % Show ( )
      end do !-- iR
      call U % Stream_SA % Show ( )
    end if !-- allocated PositionSpace_SA

  end subroutine ShowDiagnostics


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iR
    real ( KDR ), dimension ( : ), allocatable :: &
      dT_RS

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )
    select type ( I )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_H_Form )
    associate &  !-- See InitializeIntegrator subroutine herein
      ( dT_1  =>  dT_Candidate ( 1 ), &
        dT_2  =>  dT_Candidate ( 2 ), &
        dT_3  =>  dT_Candidate ( 3 ) )

    !-- Gravity step

    call U % Compute_dT_G_CGS ( dT_1, iC, T_Option )
    dT_1  =  U % GravityFactor  *  dT_1    

    !-- Fluid advection step

    if ( U % Coarsen ) then
      call U % Compute_dT_CS_CGS_C &
             ( I % EigenspeedSet_X, dT_2, U % Coarsening_F, iC, T_Option )
    else !-- .not. Coarsen
      call I % Compute_dT_CS_CGS &
             ( I % EigenspeedSet_X, dT_2, iC, T_Option )
    end if !-- Coarsen
    dT_2  =  I % CourantFactor  *  dT_2
    
    !-- Radiation streaming step

    associate ( nR  =>  U % nRadiations )
    allocate ( dT_RS ( nR ) )
    dT_RS  =  dT_3
    do iR  =  1, nR
      if ( U % Coarsen ) then
        call U % Compute_dT_CS_CGS_C &
               ( I % EigenspeedSet_X_1D ( :, iR ), dT_RS ( iR ), &
                 U % Coarsening_R ( iR ), iC, T_Option )
      else !-- .not. Coarsen
        call I % Compute_dT_CS_CGS &
               ( I % EigenspeedSet_X_1D ( :, iR ), dT_RS ( iR ), iC, T_Option )
      end if
    end do !-- iR
    end associate !-- nR
    dT_3  =  I % CourantFactor_1D  *  minval ( dT_RS )

    end associate !-- dT_1, etc.
    end select !-- S
    end select !-- I
    end select !-- U

  end subroutine Compute_dT_Local


  subroutine InitializeSeries ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    allocate ( Series_R_CC_C_Form :: I % Series )

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )
    select type ( I )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( S  =>  I % Series )
      class is ( Series_R_CC_C_Form )
    select type ( M  =>  U % Measures )
      class is ( Measures_R_CC_C_Form )
    call S % Initialize &
      ( M, I % CurrentSet_X_1D, I % CurrentSet_X, &
        I % GridImageStream, I % dT_Label, I % Unit_T, I % dT_Candidate, &
        I % T, I % Communicator % Rank, I % nWrite, I % iCycle )
    end select !-- M
    end select !-- S
    end select !-- I
    end select !-- U

  end subroutine InitializeSeries


  subroutine Analyze ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iR

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call Average ( U )

    !-- Measures
    associate ( M  =>  U % Measures )
    call M % Compute ( )
    end associate !-- M

    end select !-- U

    select type ( I )
      class is ( Integrator_CS_Form )
    call I % Analyze_CS ( I, Ignorability, T_Option )
    end select !-- I

  end subroutine Analyze


  subroutine Write ( I, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call U % Write_F_C ( I, T_Option )

    end select !-- U

  end subroutine Write


  subroutine Average ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    !-- Fluid and gravity
    call U % Average_F_C ( )

    !-- Azimuthal average
    if ( allocated ( U % PositionSpace_AA ) ) then
      do iR  =  1,  U % nRadiations
        call U % AA_Radiation ( iR ) % Compute ( )
        select type ( R_AA  =>  U % AA_Radiation ( iR ) % FieldSet_AA )
          class is ( NeutrinoMoments_G_Form )
        select type ( I_AA  =>  U % AA_Interactions ( iR ) % FieldSet_AA )
          class is ( Interactions_NM_G_Form )
        call R_AA % ComputeFromBalanced ( )
        call I_AA % Compute ( )
        end select !-- I_AA
        end select !-- R_AA
      end do !-- iR
    end if !-- allocated PositionSpace_AA

    !-- Spherical average
    if ( allocated ( U % PositionSpace_SA ) ) then
      do iR  =  1,  U % nRadiations
        call U % SA_Radiation ( iR ) % Compute ( )
        select type ( R_SA  =>  U % SA_Radiation ( iR ) % FieldSet_SA )
          class is ( NeutrinoMoments_G_Form )
        select type ( I_SA  =>  U % SA_Interactions ( iR ) % FieldSet_SA )
          class is ( Interactions_NM_G_Form )
        call R_SA % ComputeFromBalanced ( )
        call I_SA % Compute ( )
        end select !-- I_SA
        end select !-- R_SA
      end do !-- iR
    end if !-- allocated PositionSpace_SA

  end subroutine Average


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call U % Set_T_CheckpointInterval_F_CC ( I )

    end select !-- U

  end subroutine Set_T_CheckpointInterval


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


end module Universe_R_CC_C__Form
